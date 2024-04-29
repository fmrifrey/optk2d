function kdata_corr = epighostcor(kspace,kdata,delay,tol)

kxd = diff(kspace(:,:,1),1); % diff the x kspace
kxd = kxd/max(abs(kxd(:))); % normalize it

msk_right = kxd > tol;
msk_left = kxd < -tol;

kdata_corr = kdata;

for shotn = 1:size(kdata,2)
    
    I_right = find(msk_right(:,shotn));
    I_left = find(msk_left(:,shotn));
    
    if ~isempty(I_right)
        kdata_corr(I_right,shotn,:) = ...
            interp1( I_right, ...
            kdata(I_right,shotn,:), ...
            I_right + delay, ...
            'linear', 0);
    end
    
    if ~isempty(I_left)
        kdata_corr(I_left,shotn,:) = ...
            interp1( I_left, ...
            kdata(I_left,shotn,:), ...
            I_left + delay, ...
            'linear', 0);
    end
    
end

end
