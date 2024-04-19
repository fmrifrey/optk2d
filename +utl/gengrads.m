function [G,nramp] = gengrads(sys,kspace)
% sys = toppe system structure
% kspace = kspace trajectory ([N x Nshots x 3] in cm^-1)

    % get gradient waveform w/o ramps
    gam = 4257.7; % Hz/G
    dt = sys.raster*1e-6; % s
    g = 1/gam*padarray(diff(kspace,1)/dt,[1,0,0],0); % G/cm
    nramp = [0,0];

    G = 0;
    
    for shotn = 1:size(g,2)

        % get initial and final kspace and gradients
        k1 = squeeze(kspace(1,shotn,:))'; % 1 = start ramp
        k2 = squeeze(kspace(end,shotn,:))'; % 2 = end ramp

        % form path start/end points
        C1 = [zeros(1,3);k1];
        C2 = [k2;zeros(1,3)];
        
        % calculate minimum time gradients
        if norm(C1(1,:) - C1(2,:),2) > 0
            [~,~,ramp1] = minTimeGradient(C1, [], 0, 0, ...
                sys.maxGrad, sys.maxSlew, sys.raster*1e-3);
            ramp1 = padarray(ramp1,[1,0],0,'post');
        else
            ramp1 = [];
        end
        
        if norm(C2(1,:) - C2(2,:),2) > 0
            [~,~,ramp2] = minTimeGradient(C2, [], 0, 0, ...
                sys.maxGrad, sys.maxSlew, sys.raster*1e-3);
            ramp2 = padarray(ramp2,[1,0],0,'pre');
        else
            ramp2 = [];
        end
        
        % update ramp sizes
        if size(ramp1,1) > nramp(1)
            G = padarray(G,[size(ramp1,1)-nramp(1),0,0],0,'pre');
            nramp(1) = size(ramp1,1);
        else
            ramp1 = padarray(ramp1,[nramp(1)-size(ramp1,1),0,0],0,'pre');
        end
        
        % update ramp2 size
        if size(ramp2,1) > nramp(2)
            G = padarray(G,[size(ramp2,1)-nramp(2),0,0],0,'post');
            nramp(2) = size(ramp2,1);
        else
            ramp2 = padarray(ramp2,[nramp(2)-size(ramp2,1),0,0],0,'post');
        end
        
        
        
        % append the ramps
        gshot = [ramp1;squeeze(g(:,shotn,:));ramp2];
        
        % append the gradients
        if shotn == 1
            G = reshape(gshot,[],1,3);
        else
            G = cat(2,G,reshape(gshot,[],1,3));
        end
        
    end

end

