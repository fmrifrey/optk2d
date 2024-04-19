function [x_star, cost, x_set] = tvrecon(klocs,kdata,varargin)
% klocs = [Nk x Nt x Nd]
% kdata = [Nk x Nc x Nt]

    % import functions
    import rec.*
    import utl.*

    % define defaults
    defaults = struct( ...
        'fov', [], ... % image fov - Ndx1 vector
        'N', [], ... % image dimensions - Ndx1 vector
        'lam', 0, ... % lagrange multiplier for TV
        'L', [], ... % Lipschitz constant
        'type', 'l1', ... % TV semi-norm type
        'niter', 100, ... % number of iterations
        'smap', [], ... % sensitivity map [N x Nc]
        'parallelize', 0, ... % option to parallelize frame-wise recons
        'show', 1 ... % show iterations of the recon as it happens
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);
    
    % import functions
    import rec.*
    import utl.*

    % convert N and fov to row vectors
    arg.N = arg.N(:)';
    arg.fov = arg.fov(:)';

    % get number of time points
    Nt = size(klocs,2);
    Nd = size(klocs,3);

    % check SENSE map
    if isempty(arg.smap) && size(kdata,3)>1
        warning('sense map is empty, compressing data to 1 coil...');
        kdata = ir_mri_coil_compress(kdata,'ncoil',1);
        arg.smap = ones([arg.N,1]);
    elseif isempty(arg.smap) && size(kdata,3)==1
        arg.smap = ones([arg.N,1]);
    end

    % create a cell array of forward NUFFT operators for each time pt
    F = cell(Nt,1);
    nufft_args = {arg.N, 6*ones(Nd,1), 2*arg.N, arg.N/2, 'table', 2^10, 'minmax:kb'};
    for n = 1:Nt
        omega = pi^2*arg.fov./arg.N.*squeeze(klocs(:,n,:));
        F{n} = Gnufft(true(arg.N),cat(2,{omega},nufft_args)); % NUFFT
        F{n} = pipe_menon_dcf(F{n}) * F{n}; % density compensation
        F{n} = Asense(F{n},arg.smap); % sensitivity encoding
    end
    
    % define the forward and adjoint operators A and At
    A = @(x) A_fwd(x,F,arg.parallelize);
    At = @(b) A_adj(b,F,arg.parallelize);

    % calculate default L
    if isempty(arg.L)
        arg.L = pwritr(A,At,arg.N);
    end

    % initialize
    P = [];
    b = kdata;
    x_new = At(b);
    Y = x_new;
    t_k_1 = 1;
    cost = zeros(arg.niter+1,1);
    x_set = cell(arg.niter+1,1);
    res = A(x_new) - b;
    cost(1) = 1/2*norm(res(:),'fro')^2 + ...
        arg.lam*tvnorm(x_new,arg.type);
    x_set{1} = x_new;

    % loop through iterations
    for i = 1:arg.niter

        % store old values
        x_old = x_new;
        t_k = t_k_1;

        % calculate the gradient
        grad = At(A(Y) - b);
        x_new = Y - grad/arg.L;

        % denoise the new image
        [x_new,P] = tvdenoise(x_new,arg.lam,P,arg.type,arg.niter);

        % check for nans
        if any(isnan(x_new(:)))
            error('nan values in x_new... try increasing L')
        end

        % display the image
        if arg.show && mod(i,5)==0
            im(abs(x_new));
            title(sprintf('iter %d, λ = %.2g',i,arg.lam))
            drawnow
        end

        % calculate t_{k+1}
        t_k_1 = (1 + sqrt(1+4*t_k^2))/2;

        % update Y
        Y = x_new + (t_k-1)/t_k_1 * (x_new - x_old);

        % calculate residual and save cost
        res = A(x_new) - b;
        cost(i+1) = 1/2*norm(res(:),'fro')^2 + ...
            arg.lam*tvnorm(x_new,arg.type);

        % save image
        x_set{i+1} = x_new;

    end

    % save estimate
    x_star = x_new;

end

function W = pipe_menon_dcf(F,itrmax)
% F = NUFFT operator
% itrmax = number of iterations

    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end

    % initialize weights to 1 (psf)
    w = ones(size(F,1),1);
    
    % loop through iterations
    for itr = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, w) ) );
        w = w ./ d;
        
    end
    
    % normalize weights
    w = w / sum(abs(w));
    W = Gdiag(w);
end