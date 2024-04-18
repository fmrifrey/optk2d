% import functions
import rec.*
import utl.*

% set parameters
N = 256;
fov = 24;

% load image
x_gt = phantom('Modified Shepp-Logan',N);

% convert to 2d
N = N*ones(1,2);
fov = fov*ones(1,2);

% create grids
[kx_grid,ky_grid] = imgrid(N./fov,N);
[x_grid,y_grid] = imgrid(fov,N);

% simulate a sensitvity map
ncoils = 4;
if ncoils > 1
    smap = mri_sensemap_sim('nx', N(1), 'ny', N(2), ...
        'dx', fov(1)/N(1), 'dy', fov(2)/N(2), ...
        'rcoil', 10, 'ncoil', ncoils, 'chat', 0);
else
    smap = ones(N);
end

% generate the sampling trajectory
load ./kspace.mat
klocs = reshape(kspace(:,:,1:2),[],2);
omega = pi*fov./N.*klocs;

% generate the fwd operator using NUFFT
nufft_args = {N, [6,6], 2*N, N/2, 'table', 2^10, 'minmax:kb'};
A = Gnufft(true(N),cat(2,{omega},nufft_args)); % NUFFT
A = Asense(A,smap);

% NUFFT the object - inverse crime k-space signal estimation
snr = 1;
kdata = A_fwd(x_gt,{A},0);
kdata = awgn(kdata,snr);

% recon
niter = 100;
tvtype = 'l1';
lams = [0 1e-2];
costs = zeros(niter+1,length(lams));
xs = cell(niter+1,1);
for i = 1:length(lams)
    fprintf('i=%d/%d\n',i,length(lams));
    [~,costs(:,i),xs{i}] = tvrecon(reshape(klocs,[],1,2), kdata, ...
        'fov', fov, 'N', N, 'smap', smap, 'niter', niter, ...
        'lam', lams(i), 'type', tvtype, 'L', 5);
end

%% Make figures

% create a figure showing k-space sampling pattern
cfigopen('tvrecon(): k-space sampling')
subplot(1,2,1)
imagesc(x_grid(:),y_grid(:),abs(x_gt));
xlabel('x (cm)');
ylabel('y (cm)');
title('ground truth')
subplot(1,2,2)
imagesc(kx_grid(:),ky_grid(:),log(abs(fftc(x_gt))+eps()))
hold on
plot(klocs(:,1),klocs(:,2),'--r.')
xlabel('kx (cm^{-1})');
ylabel('ky (cm^{-1})');
hold off
title('k-space sampling');

% create a figure showing grid of lambdas and iterations
cfigopen('tvrecon(): λ sweep')
nRMSE = @(x) 100*sqrt(mean((x/norm(x(:)) - x_gt/norm(x_gt(:))./(x/norm(x(:)))).^2,'all'));
for i = 1:length(lams)
    x_set = xs{i};
    subplot(length(lams),3,(i-1)*3 + 1), imagesc(abs(x_set{1}))
    title(sprintf('λ = %.2g, 0 itr\nnRMSE = %.2g%%',lams(i),nRMSE(x_set{1}))); axis off
    subplot(length(lams),3,(i-1)*3 + 2), imagesc(abs(x_set{round((niter+1)/2)}))
    title(sprintf('λ = %.2g, %d itr\nnRMSE = %.2g%%',lams(i),round((niter+1)/2),nRMSE(x_set{round((niter+1)/2)}))); axis off
    subplot(length(lams),3,(i-1)*3 + 3), imagesc(abs(x_set{niter+1}))
    title(sprintf('λ = %.2g, %d itr\nnRMSE = %.2g%%',lams(i),niter+1,nRMSE(x_set{niter+1}))); axis off
end

% create a figure plotting the cost functions for each lamba
cfigopen('tvrecon(): cost functions')
plot(costs);
ylim([0 1.5*costs(1,1)])
xlabel('iteration #');
ylabel('cost');
labels = cell(length(lams),1);
for i = 1:length(lams), labels{i} = sprintf('λ = %.2g',lams(i)); end
legend(labels)