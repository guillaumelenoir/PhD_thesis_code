function [b_p_ap,b_p_am,b_m_ap,b_m_am,ap,am]=isocontour_morlet_exact(w0,R)

    % Iso-contour of the reproducing kernel at level R - EXACT solution

    theta=linspace(0,pi/2,2500);
    r=sqrt(-2*log(R));
    beta=lambertw(exp(-r^2.*sin(theta).^2+w0^2).*w0^2);         % Issue: exponential numerically too big for large w0. 
    ap=w0^2./beta+sqrt(w0^4./beta.^2-1);
    am=w0^2./beta-sqrt(w0^4./beta.^2-1);
    b_p_ap=sqrt(r^2.*cos(theta).^2.*(ap.^2+1));
    b_p_am=sqrt(r^2.*cos(theta).^2.*(am.^2+1));
    b_m_ap=-b_p_ap;
    b_m_am=-b_p_am;
    
end