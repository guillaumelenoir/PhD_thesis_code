function bound=bound_support_frequency_paul(eps)

    hat_psi=@(x) 4.*x.^2.*exp(-2.*x);    % its norm from 0 to inf = 1
    bound=1;
    while 1==1
        norm_approx=integral(hat_psi,0,bound);
        if (1-norm_approx)<eps
            break
        else
            bound=bound+1;
        end
    end
    
end