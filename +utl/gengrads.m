function [G,nramp] = gengrads(sys,kspace,rtype)
% sys = toppe system structure
% kspace = kspace trajectory ([N x Nshots x 3] in cm^-1)
% rtype = ramp type ('mint' for min time, or 'trap' for trapezoids)

    if nargin<3 || isempty(rtype)
        rtype = 'trap';
    end
        
    % get gradient waveform w/o ramps
    gam = sys.gamma*1e-4; % Hz/G
    dt = sys.raster*1e-6; % s
    g = 1/gam*padarray(diff(kspace,1)/dt,[1,0,0],0,'pre'); % G/cm
    nramp = [0,0];

    G = 0;
    
    for shotn = 1:size(g,2)

        % get initial and final kspace and gradients
        k1 = squeeze(kspace(1,shotn,:))'; % 1 = start ramp
        k2 = squeeze(kspace(end,shotn,:))'; % 2 = end ramp

        % form path start/end points
        C1 = [zeros(1,3);k1];
        C2 = [k2;zeros(1,3)];
        
        % calculate ramp 1
        if norm(k1,2) > 0 && strcmpi(rtype,'mint')
            [~,~,ramp1] = minTimeGradient(C1, [], 0, 0, ...
                sys.maxGrad, sys.maxSlew, sys.raster*1e-3);
            ramp1 = padarray(ramp1,[1,0],0,'post'); % add an extra dead sample
        elseif norm(k1,2) > 0 && strcmpi(rtype,'trap')
            ramp1 = k1/norm(k1,2).*toppe.utils.trapwave2(norm(k1,2)/gam,sys.maxGrad,sys.maxSlew,dt*1e3)';
            ramp1 = padarray(ramp1,[1,0],0,'post');
        else
            ramp1 = [];
        end
        
        % calculate ramp 2
        gf = squeeze(g(end,shotn,:));
        gf = gf(:)';
        if norm(k2,2) > 0 && strcmpi(rtype,'mint')
            [~,~,ramp2] = minTimeGradient(C2, [], gf, 0, ...
                sys.maxGrad, sys.maxSlew, sys.raster*1e-3);
            ramp2 = padarray(ramp2,[1,0],0,'pre');
        elseif norm(k2,2) > 0 && strcmpi(rtype,'trap')
            
            % ramp down the gradients
            if norm(gf,2) > 0
                ngramp = ceil(norm(gf,2)/(sys.maxSlew*1e3)/dt);
                
                % update kspace
                k2 = k2 + 0.5*dt*(ngramp-1).*gam*gf;
                
                ramp2_pre = gf.*linspace(1-1/ngramp,0,ngramp)';
            else
                ramp2_pre = [];
            end
            
            % rewind kspace with a trapezoid
            ramp2_trap = -k2/norm(k2,2).*toppe.utils.trapwave2(norm(k2,2)/gam,sys.maxGrad,sys.maxSlew,dt*1e3)';
            ramp2_trap = padarray(ramp2_trap,[1,0],0,'pre');
            
            % combine
            ramp2 = cat(1,ramp2_pre,ramp2_trap);
            
        else
            ramp2 = [];
        end
        
        % update ramp1 size
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

