function v=paul_analytical_4_eigenvectors_frequency(omega)

    % First four normalized eigenvectors

    n_omega=numel(omega);
    v=zeros(n_omega,4);
    v(:,1)=omega.*exp(-omega);
    v(:,1)=v(:,1)./sqrt(2);
    v(:,2)=omega.*exp(-omega).*(-2.*omega+2+1);
    v(:,2)=v(:,2)./sqrt(6);
    v(:,3)=omega.*exp(-omega).*0.5.*(4.*omega.^2-2.*(2+2).*2.*omega+(2+1).*(2+2));
    v(:,3)=v(:,3)./sqrt(12);
    v(:,4)=omega.*exp(-omega).*(-1/6).*(8.*omega.^3-3.*(2+3).*4.*omega.^2+3.*(2+2).*(2+3).*2.*omega-(2+1).*(2+2).*(2+3));
    v(:,4)=v(:,4)./sqrt(20);  
    v=v.*2.*sqrt(2);

end