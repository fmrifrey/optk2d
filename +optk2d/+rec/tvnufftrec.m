function [x_star, cost, x_set] = tvnufftrec(klocs,kdata,N,fov,varargin)
% klocs = kspace sampling locations [Nk x Nt x Nd]
% kdata = kspace data at sampling locs [Nk x Nc x Nt]
% N = image dimensions [Nd x 1]
% fov = field of view (inverse units of klocs), [Nd x 1]

    % define defaults
    defaults = struct( ...
        'lam', 0, ... % lagrange multiplier for TV
        'L', [], ... % Lipschitz constant
        'type', 'l1', ... % TV semi-norm type
        'niter', 100, ... % number of iterations
        'smap', [], ... % sensitivity map [N x Nc]
        'parallelize', 0, ... % option to parallelize frame-wise recons
        'denscomp', 1, ... % option to perform sampling density compensation
        'show', 0 ... % show iterations of the recon as it happens
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);

    % convert N and fov to row vectors
    N = N(:)';
    fov = fov(:)';

    % get number of time points
    Nt = size(klocs,2);
    Nd = size(klocs,3);

    % check SENSE map
    if isempty(arg.smap) && size(kdata,2)>1
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
        if arg.denscomp
            W = Gdiag(vecnorm(omega,2,2)+0.3);% optk2d.rec.pmdcf(F{n});
            F{n} = W * F{n}; % density compensation
        end
        F{n} = Asense(F{n},arg.smap); % sensitivity encoding
    end
    
    % define the forward and adjoint operators A and At
    A = @(x) optk2d.rec.A_fwd(x,F,arg.parallelize);
    At = @(b) optk2d.rec.A_adj(b,F,arg.parallelize);

    % calculate default L
    if isempty(arg.L)
        arg.L = optk2d.rec.pwritr(A,At,N);
    end

    % recon the data
    [x_star,cost,x_set] = optk2d.rec.tvrecon(A,At,kdata, ...
        'lam', arg.lam, ...
        'L', arg.L, ...
        'type', arg.type, ...
        'niter', arg.niter, ...
        'show', arg.show);

end