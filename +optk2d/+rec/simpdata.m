function [kspace,kdata,smap,x_gt] = simpdata(varargin)

    % define defaults
    defaults = struct( ...
        'kspace', [], ... % kspace sampling locations (empty to read traj from ./kspace.mat)
        'N', 256*ones(1,2), ... % isotropic matrix size
        'fov', 24*ones(1,2), ... % isotropic fov
        'ncoils', 4, ... % number of coils to simulate
        'snr', inf, ... % signal to noise ratio (see agwn.m)
        'show', 1 ... % show the ground truth image and sampling scheme
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);

    arg.N = arg.N(:)';
    arg.fov = arg.fov(:)';
    nd = size(arg.N,2);
    
    % generate phantom image & sensitivity maps
    switch nd
        case 2
            x_gt = phantom('modified shepp-logan',arg.N(1));
            smap = ir_mri_sensemap_sim('nx', arg.N(1), 'ny', arg.N(2),  ...
                'dx', arg.fov(1)/arg.N(1), 'dy', arg.fov(2)/arg.N(2), ...
                'rcoil', 10, 'ncoil', arg.ncoils, 'chat', 0);
        case 3
            x_gt = optk2d.utl.phantom3d('modified shepp-logan', arg.N(1));
            x_gt = reshape(x_gt,arg.N);
            smap = ir_mri_sensemap_sim('nx', arg.N(1), 'ny', arg.N(2), 'nz', arg.N(3), ...
                'dx', arg.fov(1)/arg.N(1), 'dy', arg.fov(2)/arg.N(2), 'dz', arg.fov(3)/arg.N(3), ...
                'rcoil', 10, 'ncoil', arg.ncoils, 'chat', 0);
        otherwise
            error('nd must be either 2 or 3 (2d or 3d)');
    end

    % load the sampling trajectory
    if isempty(arg.kspace) && isfile('./kspace.mat')
        load ./kspace.mat kspace
    elseif isempty(arg.kspace)
        error('no kspace.mat file in working dir, must specify trajectory')
    else
        kspace = arg.kspace;
    end
    
    % reformat for NUFFT
    klocs = reshape(kspace(:,:,1:nd),[],nd);
    omega = 2*pi*arg.fov./arg.N.*klocs;
    
    % generate the fwd operator using NUFFT
    nufft_args = {arg.N, 6*ones(1,nd), 2*arg.N, arg.N/2, 'table', 2^10, 'minmax:kb'};
    A = Gnufft(true(arg.N),cat(2,{omega},nufft_args)); % NUFFT
    A = Asense(A,smap);
    
    % NUFFT the object - inverse crime k-space signal estimation
    kdata = optk2d.rec.A_fwd(x_gt,{A},0);
    kdata = awgn(kdata,arg.snr);
    
    % create a figure showing k-space sampling pattern
    if arg.show && nd < 3
        optk2d.utl.cfigopen('simpdata(): k-space sampling')
        
        subplot(1,3,1)
        [x_grid,y_grid] = optk2d.utl.imgrid(arg.fov,arg.N);
        imagesc(x_grid(:),y_grid(:),abs(x_gt));
        xlabel('x (cm)');
        ylabel('y (cm)');
        title('ground truth')
        
        subplot(1,3,2)
        for i = 1:size(kspace,2)
            kdata2 = reshape(kdata,size(kspace,1),size(kspace,2),[]);
            col = mean(abs(kdata2(:,i,:)),3)./max(abs(kdata2(:,i,:)),[],'all');
            x = kspace(:,i,1)';
            y = kspace(:,i,2)';
            z = zeros(1,size(kspace,1));
            c = col';
            surf([x;x],[y;y],[z;z],[c;c], ...
                'facecol', 'no', ...
                'edgecol', 'interp', ...
                'linew', 2 ...
                );
            hold on
        end
        hold off
        view(2)
        set(gca,'color','k');
        xlabel('kx (cm^{-1})');
        ylabel('ky (cm^{-1})');
        title('k-space sampling');
        
        subplot(1,3,3)
        imagesc(x_grid(:),y_grid(:),abs(A'*kdata));
        title("A'y reconstruction");
        xlabel('x (cm)');
        ylabel('y (cm)');
        
        drawnow
    elseif arg.show
        warning('3d figures for gendata are not yet supported, sorry!\n');
    end

end

