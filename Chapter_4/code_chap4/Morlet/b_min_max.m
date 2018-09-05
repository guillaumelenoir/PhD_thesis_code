function [b_min,b_max,a_b_sq_max]=b_min_max(w0,R,first_guess)

    % Improvement required: first_guess should be done internally, because the interval to give
    % must always contains 1. 

    b_sq=@(s) -w0^2*(s-1).^2+(s.^2+1).*log(2.*s./R.^2./(s.^2+1));
    b_sq_der=@(s) s.*(-2*w0^2-1+2.*log(2.*s./R^2./(s.^2+1)))+2*w0^2+1./s;
    a_b_sq_max=fzero(b_sq_der,first_guess);
    b_max=sqrt(b_sq(a_b_sq_max));    
    b_min=-b_max;

end