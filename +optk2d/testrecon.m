% set parameters
N = 256*ones(1,2);
fov = 24*ones(1,2);
niter = 100;
lams = [0,0.05,0.1,0.5,1];
tvtype = 'l1';

% simulate data
[kspace,kdata,smap,x_gt] = optk2d.rec.simpdata('kspace', [], ... % kspace sampling locations (empty to read traj from ./kspace.mat)
    'N', N, ... % isotropic matrix size
    'fov', fov, ... % isotropic fov
    'ncoils', 32, ... % number of coils to simulate
    'snr', 1, ... % signal to noise ratio (see agwn.m)
    'show', 1 ... % show the ground truth image and sampling scheme
    );

% loop through settings
recons = cell(length(lams),1);
for i = 1:length(lams)
    
    % recon
    [x_star,cost,x_set] = optk2d.rec.tvnufftrec(reshape(kspace(:,:,1:length(fov)),[],1,length(fov)), kdata, ...
        N, fov, ...
        'smap', smap, 'niter', niter, ...
        'lam', lams(i), 'type', tvtype, 'show', 1);
    
    % save recon with info
    recons{i} = struct( ...
        'fov', fov, 'N', N, ...
        'niter', niter, 'lam', lams(i), ...
        'x_star', x_star, 'cost', cost, 'x_set', squeeze(cat(4,x_set{:})));
    
end
save test.mat recons

%% figures
nRMSE = @(x) 100*sqrt(mean((x/norm(x(:)) - x_gt/norm(x_gt(:))./(x/norm(x(:)))).^2,'all'));

% create grids
[kx_grid,ky_grid] = optk2d.utl.imgrid(N./fov,N);
[x_grid,y_grid] = optk2d.utl.imgrid(fov,N);

% create a figure showing grid of lambdas and iterations
optk2d.utl.cfigopen('testrecon(): λ sweep')
costs = zeros(niter+1,length(lams));
for i = 1:length(lams)
    x_set = recons{i}.x_set;
    
    subplot(length(lams),3,(i-1)*3 + 1)
    imagesc(abs(x_set(:,:,1)))
    title(sprintf('λ = %.2g, 0 itr\nnRMSE = %.2g%%',lams(i),nRMSE(x_set(:,:,1))));
    axis off
    
    subplot(length(lams),3,(i-1)*3 + 2)
    imagesc(abs(x_set(:,:,round((niter+1)/2))))
    title(sprintf('λ = %.2g, %d itr\nnRMSE = %.2g%%',lams(i),round((niter+1)/2),nRMSE(x_set(:,:,round((niter+1)/2))))); axis off
    
    subplot(length(lams),3,(i-1)*3 + 3)
    imagesc(abs(x_set(:,:,niter+1)))
    title(sprintf('λ = %.2g, %d itr\nnRMSE = %.2g%%',lams(i),niter+1,nRMSE(x_set(:,:,niter+1)))); axis off
    
    costs(:,i) = recons{i}.cost(:);
end

% create a figure plotting the cost functions for each lamba
optk2d.utl.cfigopen('tvrecon(): cost functions')
surf(lams,1:niter+1,costs);
zlim([min(costs,[],'all'),costs(1,1)*1.2]);
xlabel('λ')
ylabel('iteration #');
zlabel('cost');