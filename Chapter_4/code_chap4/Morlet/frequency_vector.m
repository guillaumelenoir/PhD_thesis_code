function [omega,domega,nomega]=frequency_vector(w0,a_min,a_max,b_max,eps)

    bound_support=bound_support_morlet(eps);
    t_bound=b_max+bound_support*a_max;
    domega=pi/t_bound;
    % Option 1: Impose 0 to be the lowest freq.
    %omega_bound=(bound_support+w0)/a_min;
    %nomega=ceil(omega_bound/domega)+1; 
    %omega=linspace(0,omega_bound,nomega);
    % Option 2 (better): Do not impose 0 as the lowest freq.
    omega_bound_sup=(w0+bound_support)/a_min;
    omega_bound_inf=(w0-bound_support)/a_max;       % TO DO: check whether it is <0
    nomega=ceil((omega_bound_sup-omega_bound_inf)/domega)+1;
    omega=linspace(omega_bound_inf,omega_bound_sup,nomega);
    %%%%%%%%%
    omega=omega';
    domega=omega(2)-omega(1);
    disp(['number of grid points in frequency = ',num2str(nomega)]);
    
end