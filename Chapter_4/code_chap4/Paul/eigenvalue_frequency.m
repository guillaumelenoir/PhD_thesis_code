% Morse wavelet
% Computes the eigenvalue from a given eigenvector, for any level of the iso-contour

function myeig=eigenvalue_frequency(v,C,omega,eps_time,eps_freq)

    format longG

    % parameters and recommended default values:
    % ------------------------------------------
    % v: Eigenvector for which we want the eigenvalue.
    % C: Radius of the iso-contour (which is a circle). 
    % Omega: the frequency grid of the vector v
    % eps_time and eps_freq: the precisions for the numerical support bounds. Take the same values as the ones used to compute v.
    
    tic
    
    % c_phi
    c_phi=2.*pi;
    
    % New frequency
    [omega_new,domega_new,~]=frequency_vector(C,eps_time,eps_freq);
    if omega_new(1)<=omega(1) && omega_new(end)>=omega(end)  % if R is smaller, the eigenvectors impose its bounds. 
        for k=1:numel(omega_new)
            if omega_new(k)>=omega(1)
                break
            end
        end
        myind_inf=k;
        for k=numel(omega_new):(-1):1
            if omega_new(k)<=omega(end)
                break
            end
        end
        myind_sup=k;
        omega_new=omega_new(myind_inf:myind_sup);
    end
    omega_new=omega_new';
    
    % v^ at the values of omega_new
    domega=omega(2)-omega(1);
    v_new=reconstruction_signal(v,domega,omega(1),omega_new);
    v_new=v_new./norm(v_new);
    
    % function to be integrated, in polar coordinates. Integration on b>=0.
    size_v_new=size(v_new);
    if size_v_new(1)==1
        v_new=transpose(v_new);
    end
    myfun_pos=@(theta,r) r.*(C+r.*sin(theta)).*abs((exp(1i.*omega_new.*r.*cos(theta)).*omega_new.*exp(-(C+r.*sin(theta)).*omega_new))*v_new)^2;    
    myfun_neg=@(theta,r) r.*(C+r.*sin(theta)).*abs((exp(-1i.*omega_new.*r.*cos(theta)).*omega_new.*exp(-(C+r.*sin(theta)).*omega_new))*v_new)^2;    
    myfun=@(theta,r) myfun_pos(theta,r)+myfun_neg(theta,r);

    % Numerical integration via nested integrals ("integral" is a built-in matlab function)
    myeig=integral(@(theta) integral(@(r) myfun(theta,r),0,sqrt(C^2-1),'ArrayValued',true), -pi/2, pi/2,'ArrayValued',true);
    disp(['Time spent for computing the eigenvalue= ', num2str(toc)]);
    myeig=myeig./c_phi.*domega_new.*4;
    
end
        