% load in raw data from pfile (toppe must be on path)
pf = dir('P*.7');
p = toppe.utils.loadpfile(pf(1).name);

% load in kspace trajectory
load kspace.mat

% get dimensions
ndat = size(p,1);
nshots = size(kspace,2);
ncoils = size(p,2);

% reshape into format desired by recon
klocs = reshape(kspace(1:ndat,:,1:2),[],1,2);
kdata = reshape(permute(p(:,:,:,:,1:nshots),[1,5,2,3,4]),[],ncoils);

[x,cost,x_set] = tvrecon(klocs,kdata,'fov',[24,24],'N',[128,128],'niter',100,'lam',0);