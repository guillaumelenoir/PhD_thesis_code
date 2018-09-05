function c_phi=c_phi(w0)

    % c_phi factor for the Morlet wavelet. Useful for the reconstruction
    % formula. 
    % Estimation with Monte Carlo methods.
    % w0>=5.5 : parameter of the Morlet wavelet
    
    if w0<5.5
        error('w0 must be >= 5.5');
    end
    
    myfun=@(x) 2.*sqrt(pi).*exp(-(x-w0).^2)./x;
    c_phi=integral(myfun,w0-5.5,w0+5.5);        % exp(-(x-w0).^2) is a Gaussian centered in w0 and with std=1
    
end
