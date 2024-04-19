%% set parms
fov = 24;
N = 256;
nd = 2;

%% get pfile and trajectory
p = dir('P*.7');
p = toppe.utils.loadpfile(p(1).name);
raw = permute(p,[1,5,2,3,4]);
raw = raw(end:-1:1,:,:);
load kspace.mat

%% generate SENSE map
xc = cell(size(raw,3),1);
for c = 1:size(raw,3)
    kmask = vecnorm(kspace,2,3) < 0.25*max(vecnorm(kspace,2,3),[],'all');
    xc{c} = rec.tvrecon(reshape(kspace(:,:,1:nd),[],1,nd), ...
        reshape(raw(:,:,c).*kmask,[],1), ...
        'N',N*ones(1,nd),'fov',fov*ones(1,nd),'niter',0,'L',0.1);
end
smap = mri_sensemap_denoise(cat(3,xc{:}),'niter',1);

%% recon the image
x = rec.tvrecon(reshape(kspace(:,:,1:nd),[],1,nd), reshape(raw,[],size(raw,3)), ...
    'N',N*ones(1,nd),'fov',fov*ones(1,nd),'niter',0,'L',0.1,'smap',smap);