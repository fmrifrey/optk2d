function P = L_adj(x)
% x = image of any dimensions

    % import functions
    import recon.*
    import tools.*

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % initialize cell array of finite diff matrices
    P = cell(nd,1);
    
    % loop through dimensions
    for d = 1:nd
        % calculate finite diff matrix for given dimension
        P{d} = circshift(x,1,d) - x;
    end

end

