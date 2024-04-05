function [G,nramp] = gengrads(sys,kspace)
% sys = toppe system structure
% kspace = kspace trajectory ([N x Nshots x 3] in cm^-1)

    % get gradient waveform w/o ramps
    gam = 4257.7; % Hz/G
    dt = sys.raster*1e-6; % s
    g = 1/gam*diff(kspace,1)/dt; % G/cm
    G = [];
    nramp = [0,0];

    for shotn = 1:size(g,2)

        % get initial and final kspace and gradients
        k0 = squeeze(kspace(1,shotn,:))';
        kf = squeeze(kspace(end,shotn,:))';
        gf = squeeze(g(end,shotn,:))';

        % calculate kspace ramp up
        if norm(k0,2) > 0
            ramp1 = k0/norm(k0,2).*toppe.utils.trapwave2(norm(k0,2)/gam,sys.maxGrad,sys.maxSlew,sys.raster*1e-3)';
        else
            ramp1 = [];
        end

        % ramp down the gradients at the end
        ngramp = ceil(norm(gf,2)/(sys.maxSlew*1e3)/dt);
        gramp = gf.*linspace(1-1/ngramp,0,ngramp)';
        
        % update kspace
        kf = kf + 0.5*dt*(ngramp-1).*gam*gf;

        % calculate kspace ramp down
        if norm(kf,2) > 0
            ramp2 = -kf/norm(kf,2).*toppe.utils.trapwave2(norm(kf,2)/gam,sys.maxGrad,sys.maxSlew,sys.raster*1e-3)';
        else
            ramp2 = [];
        end
        ramp2 = [gramp;ramp2];

        gshot = [ramp1;squeeze(g(:,1,:));ramp2];

        % append to the end of the array
        if ~isempty(G) && (size(gshot,1)>size(G,1))
            G = cat(1,G,zeros(size(gshot,1)-size(G,1),size(G,2),size(G,3)));
        elseif (size(gshot,1)<size(G,1))
            gshot = cat(1,gshot,zeros(size(G,1)-size(gshot,1),3));
        end
        G = cat(2,G,permute(gshot,[1,3,2]));

        % update ramp sizes
        if size(ramp1,1) > nramp(1)
            nramp(1) = size(ramp1,1);
        end
        if size(ramp2,1) > nramp(2)
            nramp(2) = size(ramp2,1);
        end

    end

end

