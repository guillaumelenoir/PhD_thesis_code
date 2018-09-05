function [a_min,a_max,a_center,a_spread]=a_min_max(w0,R,a_b_sq_max)

    b_sq=@(s) -w0^2*(s-1).^2+(s.^2+1).*log(2.*s./R.^2./(s.^2+1));    
    first_guess_2p=[a_b_sq_max a_b_sq_max*100];        % Must be improved because fzero sometimes fails to find the zero. 
    a_max=fzero(b_sq,first_guess_2p);
    a_min=1/a_max;
    a_center=(a_min+a_max)/2;
    a_spread=a_max-a_min;
    
end
    
