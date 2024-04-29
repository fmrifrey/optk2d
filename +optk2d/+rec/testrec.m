function [outputArg1,outputArg2] = recon(varargin)

    % load in kspace trajectory
    load kspace.mat

    % load in pfile
    p = dir('P*.7');
    if isempty(p)
        fprintf('no pfile found in %s\n--> simulating phantom data\n',pwd);
    else
        fprintf('loading 
        p = toppe.utils.loadpfile(p(1).name);
        raw = permute(p,[1,5,2,3,4]);
        raw = raw(end:-1:1,:,:);
    end

end

