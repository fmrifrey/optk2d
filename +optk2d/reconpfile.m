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

% correct sizes
if size(raw,2) > size(kspace,2)
    raw = raw(:,2:end,:);
end

if size(kspace,1) > size(raw,1)
    warning('size(kspace,1) < size(raw,1)');
    kspace = kspace(1:size(raw,1),:,:);
elseif size(kspace,1) < size(raw,1)
    warning('size(kspace,1) < size(raw,1)');
    raw = raw(1:size(kspace,1),:,:);
end

epicorr = 1;
if epicorr
    del = -0.95;
    raw = optk2d.utl.epighostcor(kspace,raw,del,1e-2);
end

%% generate SENSE map
if ~isvar('smap')
    xc = cell(size(raw,3),1);
    for c = 1:size(raw,3)
        kmask = vecnorm(kspace,2,3) < 0.25*max(vecnorm(kspace,2,3),[],'all');
        xc{c} = optk2d.rec.tvnufftrec(reshape(kspace(:,:,1:nd),[],1,nd), ...
            reshape(raw(:,:,c).*kmask,[],1), ...
            N*ones(1,nd), ...
            fov*ones(1,nd), ...
            'niter',0,'L',0.1);
    end
    smap = bart('ecalib -b0 -m1', optk2d.utl.fftc(cat(4,xc{:}),1:3));
    smap = squeeze(smap);
end
save smap.mat

%% recon the image
niter = 500;
tvtype = 'l1';
lam = 0.1;

x_nonsense = optk2d.rec.tvnufftrec(reshape(kspace(:,:,1:nd),[],1,nd), reshape(raw,[],size(raw,3)), ...
    N*ones(1,nd), ...
    fov*ones(1,nd), ...
    'smap', [], ...
    'niter', 0, ...
    'lam', 0, ...
    'L', 1, ...
    'show', 1);

% recon
[x_star,cost,x_set] = optk2d.rec.tvnufftrec(reshape(kspace(:,:,1:nd),[],1,nd), reshape(raw,[],size(raw,3)), ...
    N*ones(1,nd), ...
    fov*ones(1,nd), ...
    'smap', smap, ...
    'niter', niter, ...
    'type', tvtype, ...
    'lam', lam, ...
    'show', 1);
save recon.mat x_nonsense x_star costs x_set