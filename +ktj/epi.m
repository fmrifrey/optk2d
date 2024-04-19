function ktraj = epi(sys,fov,N,varargin)
% sys = toppe system structure
% fov = field of view (cm, scalar - assuming isotropic 2D)
% N = matrix size (scalar - assuming isotropic 2D)

    % define defaults
    defaults = struct( ...
        'Ns', [], ... % number of shots
        'ya', 1, ... % y acceleration (downsampling) factor
        'save', (nargout < 1) ... % saves an .h5 file in the current directory
        );

    % parse arguments
    arg = vararg_pair(defaults,varargin);
    
    % set number of phase encodes
    Np = ceil(N/arg.ya);
    
    % default - non epi (nshots = nphases)
    if isempty(arg.Ns)
        arg.Ns = Np;
    elseif mod(Np,arg.Ns)
        warning('Np % Ns > 0 -> planes may not be fully sampled\n');
    end
    
    % generate cartesian x end points for each freq. encode
    [Cx,Cy] = utl.imgrid(N/fov,[N,Np]);
    
    % calculate path points for first shot
    Cx_shot0 = Cx(:,1:arg.Ns:end); % index every Ns f encode
    Cy_shot0 = Cy(:,1:arg.Ns:end);
    Cx_shot0(:,1:2:end) = Cx_shot0(end:-1:1,1:2:end); % flip every other f encode
    Cy_shot0(:,1:2:end) = Cy_shot0(end:-1:1,1:2:end);
    C_shot0 = [Cx_shot0(:),Cy_shot0(:),zeros(size(Cx_shot0(:)))];
    
    % determine grad limit based on fov
    Gmax_fov = 1/(sys.gamma*1e-7) * 1/fov/(sys.raster*1e-3);
    kshot0 = minTimeGradient(C_shot0, [], 0, 0, ...
        min(Gmax_fov,sys.maxGrad), sys.maxSlew, sys.raster*1e-3);
    
    % initialize ktraj
    ktraj = zeros(size(kshot0,1),arg.Ns,3);
    
    % separate into shots
    for shotn = 1:arg.Ns
        kshotn = kshot0;
        kshotn(:,2) = kshotn(:,2) + N/Np/fov * floor((shotn-1)/2); % add shot offset
        if mod(shotn,2) % odd shots - flip x and y
            kshotn = -kshotn;
        end
        ktraj(:,shotn,:) = reshape(kshotn,[],1,3);
    end
    
    % save trajectory
    if arg.save
        h5create('./traj.h5','/ktraj',size(ktraj));
        h5write('./traj.h5','/ktraj',ktraj);
    end
    
end

