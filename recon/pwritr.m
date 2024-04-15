function L = pwritr(A,At,sz,tol)
% A = forward operator @(x) A*x
% At = adjoint operator @(y) A'*y
% sz = size of input image [Nd x 1]
% tol = tolerance

    if nargin < 4 || isempty(tol)
        tol = 1e-2;
    end

    % generate a random image
    b_k = randn(sz);
    
    for i = 1:10

        % calculate A'A*b_k
        b_k1 = At(A(b_k));
        
        % calculate the spectral radius
        L = real(b_k'*b_k1/(b_k'*b_k));
        
        % renormalize the vector
        b_k = b_k1 / norm(b_k1);

    end

end

