function [x_el,y_el]=ellipse_a_b(b_max,a_center,a_spread)    

    axe_horiz_ellipse=b_max;
    axe_vert_ellipse=a_spread/2;
    z=linspace(0,2*pi,10000);
    x_el=axe_horiz_ellipse.*cos(z);
    y_el=a_center+axe_vert_ellipse.*sin(z);
    
end