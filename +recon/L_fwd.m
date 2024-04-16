function x = L_fwd(P)
% P = cell array of projection matrices for each dimension

    % import functions
    import recon.*
    import tools.*
    
    % get size
    sz = size(P{1});
    nd = ndims(P);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % initialize image
    x = zeros(sz);

    % loop through dimensions
    for d = 1:nd

        % calculate L(p,q)
        % L(p,q)_i,j = p_i,j + q_i,j - p_i-1,j - q_i-1,j
        x = x + circshift(P{d},d) - P{d};
        
    end

end