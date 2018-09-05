function bound=bound_support_time_paul(eps)

    abs_psi=@(x) 2./pi.*(1+x.^2).^(-2); % its norm from 0 to inf = 1
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