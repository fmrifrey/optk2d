function b = A_fwd(x,F,parallelize)
% x = image [N x N x Nt]
% F = NUFFT operators {[Nk x N^2] x Nt}
% parallelize = option to parallelize the frame loop (0 or 1)

    if nargin < 3 | isempty(parallelize)
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

    b = zeros(Nk,Nc,Nt);
    parfor (n = 1:Nt, 100*parallelize) % loop through time points
        b(:,:,n) = F{n}*x(:,:,n); % calculate NUFFT of image at time point n
    end

end