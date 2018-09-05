function [omega,domega,nomega]=frequency_vector(C,eps_time,eps_freq)

    radius=sqrt(C^2-1);
    b_max=radius;
    a_min=C-radius;
    a_max=C+radius;
    t_bound=b_max+bound_support_time_paul(eps_time)*a_max;
    domega=pi/t_bound;
    omega_bound=bound_support_frequency_paul(eps_freq)/a_min;
    nomega=ceil(omega_bound/domega)+1; 
    omega=linspace(0,omega_bound,nomega);
    omega=omega';
    domega=omega(2)-omega(1);
    disp(['number of grid points in frequency = ',num2str(nomega)]);
    
end