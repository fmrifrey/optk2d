function x = L_fwd(P)
% P = cell array of projection matrices for each dimension

    % import functions
    import optk2d.rec.*
    import optk2d.utl.*
    
    % get size
    sz = size(P{1});
    sz(1) = sz(1)+1;
    nd = ndims(P);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % initialize image
    x = zeros(sz);
    idcs1 = cell(nd,1);
    idcs2 = cell(nd,1);

    % loop through dimensions
    for d1 = 1:nd
        
        % get indicies for addition and subtraction
        for d2 = 1:nd
            if d1 == d2
                idcs1{d2} = 1:sz(d2)-1;
                idcs2{d2} = 2:sz(d2);
            else
                idcs1{d2} = 1:sz(d2);
                idcs2{d2} = 1:sz(d2);
            end
        end
        
        % calculate
        % L(p,q)_i,j = p_i,j + q_i,j - p_i-1,j - q_i-1,j
        x(idcs1{:}) = x(idcs1{:}) + P{d1};
        x(idcs2{:}) = x(idcs2{:}) - P{d1};
        
    end

end