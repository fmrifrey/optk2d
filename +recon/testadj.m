function ratio = testadj(A_fwd,A_adj,xtmp)
% A_fwd = forward operator to test
% A_adj = adjoint operator to test
% xtmp = template in shape of x

if nargin<3 || isempty(xtmp)
    xtmp = zeros(64);
end

% generate random x and compute Ax
if iscell(xtmp) % for cases that x is a cell
    x = cell(size(xtmp));
    for i = 1:length(x)
        x{i} = randn(size(xtmp{i}));
    end
else
    x = randn(size(xtmp));
end
Ax = A_fwd(x);

% generate random y and compute A'y
if iscell(Ax) % for cases that y is a cell
    y = cell(size(Ax));
    for i = 1:length(y)
        y{i} = randn(size(Ax{i}));
    end
else
    y = randn(size(Ax));
end
Aty = A_adj(y);

% calculate numerator <Ax,y>
if iscell(y) % for cases that y is a cell
    numer = 0;
    for i = 1:length(y)
        numer = numer + Ax{i}(:)'*y{i}(:);
    end
else
    numer = Ax(:)'*y(:);
end

% calculate denominator <x,A'y>
if iscell(x) % for cases that y is a cell
    denom = 0;
    for i = 1:length(x)
        denom = denom + x{i}(:)'*Aty{i}(:);
    end
else
    denom = x(:)'*Aty(:);
end

% calculate <Ax,y>/<x,A'y>
ratio = numer/denom;

end

