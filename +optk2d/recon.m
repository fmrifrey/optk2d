function recon(varargin)

    % load in kspace trajectory
    load kspace.mat

    % load in pfile
    p = dir('P*.7');
    if isempty(p)
        fprintf('recon(): no pfile found in %s\n--> simulating phantom data\n',pwd);
        kdata = optk2d.rec.simpdata('kspace',kspace,'N',arg.N,'fov', ...
            arg.fov,'nd',arg.nd,'ncoils',arg.ncsim);
    else
        fprintf('recon(): loading data from pfile: %s...\n', p(1).name);
        p = toppe.utils.loadpfile(p(1).name);
        kdata = permute(p,[1,5,2,3,4]);
        kdata = kdata(end:-1:1,:,:);
    end

end

