% Morlet wavelet
% Eigenvalues/vectors of the integrated projection inside an Iso-contour of the norm of the reproducing kernel

function [omega,d,V]=integrated_projection_simple_integral_exact_contour_frequency(w0,R,eps,ind_fig)

    format longG

    % parameters and recommended default values:
    % ------------------------------------------
    % w0=5.5: parameter of the Morlet wavelet. w0 must be >= 5.5.
    % R=0.5 is the level of the contour. 0<R<1. 
    % eps=10^(-6): Value of the wavelet which is used to define the max. of the support of a Gaussian (in time or frequency)
    % ind_fig: index for the figure of the contour drawn whitin this code.
    
    % Extrema of the contour at level R
    [b_min,b_max,a_b_sq_max]=b_min_max(w0,R,[0.000001 1000000]);
    [a_min,a_max,~,~]=a_min_max(w0,R,a_b_sq_max);  
    %b_min=0;
    %r=sqrt(-2*log(R));
    %beta=lambertw(exp(-r^2+w0^2+log(w0^2)));
    %a_min=w0^2/beta-sqrt(w0^4/beta^2-1);
    %a_max=1/a_min;
    
    % Plot the contour and its extrema
    [b_p_ap,b_p_am,b_m_ap,b_m_am,ap,am]=isocontour_morlet_exact(w0,R);
    figure(ind_fig);
    hold on
    plot(b_p_ap,ap,'r','linewidth',3)
    plot(b_m_ap,ap,'b','linewidth',3)
    plot(b_p_am,am,'c','linewidth',3)
    plot(b_m_am,am,'g','linewidth',3)
    plot(0,a_min,'k.','markersize',20)
    plot(0,a_max,'k.','markersize',20)
    plot(b_max,a_b_sq_max,'k.','markersize',20)
    plot(b_min,a_b_sq_max,'k.','markersize',20)
    
    % frequency vector
    [omega,domega,nomega]=frequency_vector(w0,a_min,a_max,b_max,eps);

    % Function to be integrated (b>=0 and a<=1)
    b_contour=@(a) sqrt(-w0^2.*(a-1).^2+(a.^2+1).*log(2.*a./R.^2./(a.^2+1)));
    myfun1=@(u,v,a) 1./a.*(sin(b_contour(a).*(v-u)).*exp(-0.5.*(u.*a-w0).^2-0.5.*(v.*a-w0).^2)+sin(b_contour(1./a).*(v-u)).*exp(-0.5.*(u./a-w0).^2-0.5.*(v./a-w0).^2));
    myfun2=@(u,a) b_contour(a)./a.*(exp(-(u.*a-w0).^2)+1./a.*exp(-(u./a-w0).^2));
    
    % Numerical integration via integral2 (built-in matlab function)
    myproj=zeros(nomega,nomega);
    tic
    for k=1:nomega;
        disp([num2str(k),'/',num2str(nomega)])
        omega_k=omega(k);
        for l=1:(k-1)
            omega_l=omega(l);
            myfun_inloop=@(a) myfun1(omega_k,omega_l,a); 
            myproj(k,l)=integral(myfun_inloop,a_min,1)./(omega_l-omega_k);
            myproj(l,k)=myproj(k,l);
        end
        myfun_inloop=@(a) myfun2(omega_k,a);
        myproj(k,k)=integral(myfun_inloop,a_min,1);
    end
    disp(['Time spent= ', num2str(toc)]);
    myproj=myproj.*2./sqrt(pi)./c_phi(w0).*domega;
    
    % Eigenvalues and eigenvectors of the integrated projection
    [V,D]=eig(myproj);
    d=diag(D);

end
        