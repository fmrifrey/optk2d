function x = A_adj(b,F,parallelize)
% b = kspace data [Nk x Nc x Nt]
% F = NUFFT operators {[Nk x N^Nd] x Nt}
% parallelize = option to parallelize the frame loop (0 or 1)
    
    % import functions
    import optk2d.rec.*
    import optk2d.utl.*

    if nargin < 3 || isempty(parallelize)
        parallelize = 0;
    end
    
    % Get sizes
    N = F{1}.idim(1);
    Nk = F{1}.odim(1);
    if length(F{1}.odim) > 1
        Nc = F{1}.odim(2);
    else
        Nc = 1;
    end
    Nt = length(F);

    x = zeros(N,N,Nt);
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        xn = F{n}'*b(:,:,n); % calculate inverse NUFFT of ksignal at time point n
        x(:,:,n) = reshape(xn,F{n}.idim);
    end

end