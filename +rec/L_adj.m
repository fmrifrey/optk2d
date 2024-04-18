function P = L_adj(x)
% x = image of any dimensions

    % import functions
    import rec.*
    import utl.*

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
        % get pad array
        padsz = zeros(1,nd);
        padsz(d) = 1;

        % calculate L'(x) = {p,q}
        % p_i,j = x_i,j - x_i+1,j
        % q_i,j = x_i,j - x_i,j+1
        P{d} = padarray(-diff(x,1,d),padsz,0,'post'); % neumann bndry cond - dx/d_end = 0
    end

end

