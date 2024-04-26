function L = pwritr(A,At,sz,tol)
% A = forward operator @(x) A*x
% At = adjoint operator @(y) A'*y
% sz = size of input image [Nd x 1]
% tol = tolerance

    % import functions
    import rec.*
    import utl.*

    if nargin < 4 || isempty(tol)
        tol = 1e-2;
    end

    b_k = randn(sz);
    Ab_k = At(A(b_k));
    b_k_1 = b_k;
    norm_b_k = norm(Ab_k,'fro');
    
    while 1
    
        b_k = Ab_k/norm_b_k;
        Ab_k = (At(A(b_k)));
        norm_b_k_1 = norm(Ab_k,'fro');
        if norm(b_k_1-b_k,'fro')/norm(b_k,'fro') <= tol
            break
        else
            norm_b_k = norm_b_k_1;
            b_k_1 = b_k;
        end
    
        
    
    end
    b = b_k;
    L = real(b(:)'*Ab_k(:)/(b(:)'*b(:)));
end

