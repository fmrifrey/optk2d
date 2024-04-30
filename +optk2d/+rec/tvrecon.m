function [x_star,cost,x_set] = tvrecon(A,At,b,varargin)
% A = forward operator @(x) A*x
% At = adjoint operator @(y) A'*y
% b = data to reconstruct

    % define defaults
    defaults = struct( ...
        'lam', 0, ... % lagrange multiplier for TV
        'L', [], ... % Lipschitz constant
        'type', 'l1', ... % TV semi-norm type
        'niter', 100, ... % number of iterations
        'show', 0 ... % show iterations of the recon as it happens
        );
    
    % parse arguments
    arg = vararg_pair(defaults,varargin);
    
    % import functions
    import optk2d.rec.*
    import optk2d.utl.*
    
    % initialize
    P = [];
    x_new = At(b);
    Y = x_new;
    t_k_1 = 1;
    cost = zeros(arg.niter+1,1);
    x_set = cell(arg.niter+1,1);
    res = A(x_new) - b;
    cost(1) = 1/2*norm(res(:),'fro')^2 + ...
        arg.lam*tvnorm(x_new,arg.type);
    x_set{1} = x_new;

    % loop through iterations
    for i = 1:arg.niter

        % store old values
        x_old = x_new;
        t_k = t_k_1;

        % calculate the gradient
        grad = At(A(Y) - b);
        x_new = Y - grad/arg.L;
        
        % denoise the new image
        if abs(arg.lam) > 0
            [x_new,P] = tvdenoise(x_new,arg.lam/arg.L,P,arg.type,arg.niter);
        end

        % check for nans
        if any(isnan(x_new(:)))
            error('nan values in x_new... try increasing L')
        end

        % display the image
        if arg.show %&& mod(i,5)==0
            cfigopen('tvrecon(): iterative reconstruction');
            subplot(2,2,1)
            im(abs(x_new));
            colormap parula
            title(sprintf('iter %d, λ = %.2g',i,arg.lam))
            subplot(2,2,2)
            im(abs(optk2d.utl.fftc(x_new)));
            colormap parula
            title(sprintf('iter %d, λ = %.2g',i,arg.lam))
            subplot(2,1,2)
            plot(cost(1:i))
            title('cost')
            drawnow
        end

        % calculate t_{k+1}
        t_k_1 = (1 + sqrt(1+4*t_k^2))/2;

        % update Y
        Y = x_new + (t_k-1)/t_k_1 * (x_new - x_old);

        % calculate residual and save cost
        res = A(x_new) - b;
        cost(i+1) = 1/2*norm(res(:),'fro')^2 + ...
            arg.lam*tvnorm(x_new,arg.type);
        
        % save image
        x_set{i+1} = x_new;

    end

    % save estimate
    x_star = x_new;
    
end