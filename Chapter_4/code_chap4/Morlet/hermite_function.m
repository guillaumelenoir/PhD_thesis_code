function psi_2=hermite_function(n,x)

    % Computes the hermite function of order n from the recursion relation
    
    psi_0=(pi^(-1/4)).*exp(-0.5.*x.^2);
    psi_1=(pi^(-1/4)).*sqrt(2).*x.*exp(-0.5.*x.^2);
    
    if n==0
        psi_2=psi_0;
    elseif n==1
        psi_2=psi_1;
    else
        for k=2:n
            psi_2=(x.*psi_1-sqrt((k-1)/2).*psi_0)./sqrt(k/2);
            psi_0=psi_1;
            psi_1=psi_2;
        end
    end
    %psi_2=psi_2./norm(psi_2);
    
end