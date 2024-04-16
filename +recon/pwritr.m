function L = pwritr(A,At,sz,tol)
% A = forward operator @(x) A*x
% At = adjoint operator @(y) A'*y
% sz = size of input image [Nd x 1]
% tol = tolerance

    % import functions
    import recon.*
    import tools.*

    if nargin < 4 || isempty(tol)
        tol = 1e-2;
    end

    b_k = randn(sz);
    Ab_k = real(At(A(b_k)));
    norm_b_k = norm(Ab_k);

    b = b_k;
    L = real(b(:)'*Ab_k(:)/(b(:)'*b(:)));
    
    while 1
    
        b_k = Ab_k/norm_b_k;
        Ab_k = real(At(A(b_k)));
        norm_b_k_1 = norm(Ab_k);
        if norm(norm_b_k_1-norm_b_k) <= tol
            break
        else
            norm_b_k = norm_b_k_1;
        end
    
        b = b_k;
        L = real(b(:)'*Ab_k(:)/(b(:)'*b(:)));
    
    end

end

