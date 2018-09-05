function bound=bound_support_morlet(eps)

    % Idem in time and frequency (the Fourier transform of a Gaussian is
    % still a Gaussian!)

    abs_psi=@(x) 1./sqrt(pi).*exp(-x.^2); % its norm from 0 to inf = 1
    bound=1;
    while 1==1
        norm_approx=integral(abs_psi,-bound,bound);
        if (1-norm_approx)<eps
            break
        else
            bound=bound+1;
        end
    end
    
end