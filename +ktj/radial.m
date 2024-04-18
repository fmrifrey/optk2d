function ktraj = radial(sys,fov,N,varargin)
% sys = toppe system structure
% fov = field of view (cm, scalar - assuming isotropic 2D)
% N = matrix size (scalar - assuming isotropic 2D)

    gam = 4257.7; % Hz/G

    % define defaults
    defaults = struct( ...
        'Ns', N, ... % number of spokes
        'dir', 1, ... % spoke direction (1 = center-out, 2 = out-center, 3 = out-out)
        'save', 1 ... % saves an .h5 file in the current directory
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);

    % define initial spoke
    switch arg.dir
        case 1 % center-out
            C = [0,0,0;
                0.5,0,0];
            dtheta = 360/arg.Ns;
        case 2 % out-center
            C = [0.5,0,0;
                0,0,0];
            dtheta = 360/arg.Ns;
        case 3 % out-out
            C = [-0.5,0,0;
                0.5,0,0];
            dtheta = 180/arg.Ns;
        otherwise
            error('invalid direction');
    end

    % determine grad limit based on fov
    Gmax_fov = 1/gam * 1/fov/(sys.raster*1e-6);

    % calculate first spoke with 
    spoke0 = minTimeGradient(N/fov*C, [], 0, 0, ...
        min(sys.maxGrad,Gmax_fov), sys.maxSlew, sys.raster*1e-3);

    % calculate trajectory for each spoke
    ktraj = zeros(size(spoke0,1),arg.Ns,3);
    for i = 1:arg.Ns
        spokei = spoke0*rotz(i*dtheta)';
        ktraj(:,i,:) = reshape(spokei,[],1,3);
    end

end

