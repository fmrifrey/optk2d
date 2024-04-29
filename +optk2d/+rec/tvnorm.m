function tv = tvnorm(x,type)
% x = image of any dimensions
% type = type of tv semi-norm ('iso' or 'l1')

    % set default type
    if nargin < 2 || isempty(type)
        type = 'l1';
    end

    switch type
        case 'iso'
            tv = tvnorm_iso(x);
        case 'l1'
            tv = tvnorm_l1(x);
        otherwise
            error('invalid type: %s',type);
    end

end

function tv = tvnorm_iso(x)

    % import functions
    import optk2d.rec.*
    import optk2d.utl.*

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % calculate projection operators
    P = L_adj(x);

    % loop through dimensions
    D = zeros(sz); % sum of squares matrix
    for d = 1:nd
        % add its square to D
        D = D + P{d}.^2;
    end

    % calculate norm
    tv = sum(sqrt(D),'all');

end

function tv = tvnorm_l1(x)

    % import functions
    import optk2d.rec.*
    import optk2d.utl.*

    % get size
    sz = size(x);
    nd = ndims(x);
    if (nd==2 && sz(2)==1)
        nd = 1;
    end

    % calculate projection operators
    P = L_adj(x);

    % loop through dimensions
    tv = 0;
    for d = 1:nd
        % add l1 norm
        tv = tv + sum(abs(P{d}),'all');
    end
    
end