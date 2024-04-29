function [W,w] = pmdcf(F,itrmax)
% pipe-menon density compensation function
% F = NUFFT operator
% itrmax = number of iterations

    if nargin < 2 || isempty(itrmax)
        itrmax = 15;
    end
    
    % initialize weights to 1 (psf)
    w = ones(size(F,1),1);
    
    % loop through iterations
    for i = 1:itrmax
        
        % Pipe algorithm: W_{i+1} = W_{i} / (G * (G' * W_{i}))
        d = real( F.arg.st.interp_table(F.arg.st, ...
            F.arg.st.interp_table_adj(F.arg.st, w) ) );
        w = w ./ d;
        w = w / norm(w,1);
        W = Gdiag(w);
        
    end
    
end