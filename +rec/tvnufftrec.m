function [x_star, cost, x_set] = tvnufftrec(klocs,kdata,N,fov,varargin)
% klocs = kspace sampling locations [Nk x Nt x Nd]
% kdata = kspace data at sampling locs [Nk x Nc x Nt]
% N = image dimensions [Nd x 1]
% fov = field of view (inverse units of klocs), [Nd x 1]

    % import functions
    import rec.*
    import utl.*

    % define defaults
    defaults = struct( ...
        'lam', 0, ... % lagrange multiplier for TV
        'L', [], ... % Lipschitz constant
        'type', 'l1', ... % TV semi-norm type
        'niter', 100, ... % number of iterations
        'smap', [], ... % sensitivity map [N x Nc]
        'parallelize', 0, ... % option to parallelize frame-wise recons
        'show', 0 ... % show iterations of the recon as it happens
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);
    
    % import functions
    import rec.*
    import utl.*

    % convert N and fov to row vectors
    N = N(:)';
    fov = fov(:)';

    % get number of time points
    Nt = size(klocs,2);
    Nd = size(klocs,3);

    % check SENSE map
    if isempty(arg.smap) && size(kdata,3)>1
        warning('sense map is empty, compressing data to 1 coil...');
        kdata = ir_mri_coil_compress(kdata,'ncoil',1);
        arg.smap = ones([N,1]);
    elseif isempty(arg.smap) && size(kdata,3)==1
        arg.smap = ones([N,1]);
    end

    % create a cell array of forward NUFFT operators for each time pt
    F = cell(Nt,1);
    nufft_args = {N, 6*ones(Nd,1), 2*N, N/2, 'table', 2^10, 'minmax:kb'};
    for n = 1:Nt
        omega = 2*pi*fov./N.*squeeze(klocs(:,n,:));
        F{n} = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
        F{n} = pipe_menon_dcf(F{n}) * F{n}; % density compensation
        F{n} = Asense(F{n},arg.smap); % sensitivity encoding
    end
    
    % define the forward and adjoint operators A and At
    A = @(x) A_fwd(x,F,arg.parallelize);
    At = @(b) A_adj(b,F,arg.parallelize);

    % calculate default L
    if isempty(arg.L)
        arg.L = 1.1*pwritr(A,At,N);
    end

    % recon the data
    [x_star,cost,x_set] = tvrecon(A,At,kdata, ...
        'lam', arg.lam, ...
        'L', arg.L, ...
        'type', arg.type, ...
        'niter', arg.niter, ...
        'show', arg.show);

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