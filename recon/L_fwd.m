function x = L_fwd(P)
% P = cell array of projection matrices for each dimension

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
        x = -(circshift(x,1,d) - P{d});
    end

end