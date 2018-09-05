% Morse wavelet
% Eigenvalues/vectors of the integrated projection inside an Iso-contour of the norm of the reproducing kernel

function [omega,d,V]=integrated_projection_simple_integral_frequency(C,eps_time,eps_freq)

    format longG

    % parameters and recommended default values:
    % ------------------------------------------
    % C=2.0: Radius of the iso-contour (which is a circle). 
    % eps_time=10^(-3): Value of the wavelet which is used to define the max. of the temporal support. 
    % eps_freq=10^(-3): Value of the Fourier transform of the wavelet which is used to define the max. of the frequency support. 
    
    % frequency vector
    [omega,domega,nomega]=frequency_vector(C,eps_time,eps_freq);
        
    % c_phi
    c_phi=2.*pi;
    
    % Function to be integrated (b>=0 and a<=1)   
    b_contour=@(a) sqrt(2.*a.*C-a.^2-1);
    myfun1=@(u,v,a) 1./a.*u.*v.*(a.^2.*exp(-a.*(u+v)).*sin(b_contour(a).*(v-u))+1./a.^2.*exp(-(u+v)./a).*sin(b_contour(1./a).*(v-u)));
    myfun2=@(u,a) 1./a.*b_contour(a).*u.^2.*(a.^2.*exp(-2.*a.*u)+1./a.^3.*exp(-2.*u./a));
    
    % plot the contour
    %figure(1)
    %mya=(C-sqrt(C^2-1)):0.0001:(1/(C-sqrt(C^2-1)));
    %plot(b_contour(mya),mya)
    
    % Numerical integration via integral2 (built-in matlab function)
    myproj=zeros(nomega,nomega);
    tic
    for k=1:nomega;
        disp([num2str(k),'/',num2str(nomega)])
        omega_k=omega(k);
        for l=1:(k-1)
            omega_l=omega(l);
            myfun_inloop=@(a) myfun1(omega_k,omega_l,a); 
            myproj(k,l)=integral(myfun_inloop,C-sqrt(C^2-1),1)./(omega_l-omega_k);
            myproj(l,k)=myproj(k,l);
        end
        myfun_inloop=@(a) myfun2(omega_k,a);
        myproj(k,k)=integral(myfun_inloop,C-sqrt(C^2-1),1);
    end
    disp(['Time spent= ', num2str(toc)]);
    myproj=myproj.*8./c_phi.*domega;
    
    % Eigenvalues and eigenvectors of the integrated projection
    [V,D]=eig(myproj);
    d=diag(D);

end
        