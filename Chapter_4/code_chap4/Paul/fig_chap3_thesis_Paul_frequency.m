% Figures for chapter 3 of PhD thesis.

%%%%%%%%%%%%%%%%%%%%
%   Paul Wavelet   %
%%%%%%%%%%%%%%%%%%%%

% Iso-contours
C=[1.1, 2, 3, 4, 5, 6];
n_C=numel(C);
theta=0:0.001:2*pi;

h=figure(1);
hold on
for k=1:n_C
    x_C=sqrt(C(k)^2-1).*cos(theta);
    y_C=C(k)+sqrt(C(k)^2-1).*sin(theta);
    plot(x_C,y_C,'k','linewidth',1)
    plot(sqrt(C(k)^2-1),C(k),'r.','markersize',20)
    plot(-sqrt(C(k)^2-1),C(k),'r.','markersize',20)
    text(sqrt(C(k)^2-1)*0.1,(C(k)+sqrt(C(k)^2-1))+0.2,strcat('C=',num2str(C(k))),'fontsize',15)
end
box off
grid on
axis([-6.5 6.5 0 13])
xlabel('b','fontsize',20)
ylabel('a','fontsize',20)
set(gca, 'FontSize', 15)
axis square     % implies that a circle will look like a circle!
saveas(h,'iso_contours_paul.eps','epsc')
close all

% 4 first eigenvectors: numerical vs analytical solutions
C=4;                              
eps_time=0.001;                    
eps_freq=0.001;                  
[omega_num,d,V]=integrated_projection_simple_integral_frequency(C,eps_time,eps_freq);
domega_num=omega_num(2)-omega_num(1);
omega=linspace(0,omega_num(end),5000);      % numel(omega_num)
omega=omega';
domega=omega(2)-omega(1);
eigenvec_num=zeros(numel(omega),4);
for k=1:4
    eigenvec_num(:,end-k+1)=reconstruction_signal(V(:,end-k+1),domega_num,0,omega);
    eigenvec_num(:,end-k+1)=eigenvec_num(:,end-k+1)./norm(eigenvec_num(:,end-k+1))./sqrt(domega);    
end
eigenvec_anal=paul_analytical_4_eigenvectors_frequency(omega);
h=figure(1);
loglog(omega,abs(eigenvec_anal(:,4)-eigenvec_num(:,end-3)),'k','linewidth',1)
hold on
loglog(omega,abs(eigenvec_anal(:,3)-eigenvec_num(:,end-2)),'m','linewidth',1)
loglog(omega,abs(eigenvec_anal(:,2)+eigenvec_num(:,end-1)),'r','linewidth',1)
loglog(omega,abs(eigenvec_anal(:,1)-eigenvec_num(:,end)),'b','linewidth',1)
lgd=legend(num2str(d(end),'%6.5g'),num2str(d(end-1),'%6.5g'),num2str(d(end-2),'%6.5g'),num2str(d(end-3),'%6.5g'));
set(lgd,'FontSize',15); 
box off
grid off
xlim([0 20])
xlabel('\omega','fontsize',25)
title('Absolute error','fontsize',20)
set(gca, 'FontSize', 15)
saveas(h,'eigenvec_paul_error.eps','epsc')
close all
h=figure(1);
semilogx(omega,eigenvec_anal(:,1),'b','linewidth',4)
hold on
semilogx(omega,eigenvec_num(:,end),'r','linewidth',2)
lgd=legend('Analytical','Numerical');
set(lgd,'FontSize',20); 
box off
grid off
xlim([0 20])
xlabel('\omega','fontsize',25)
title(strcat('\lambda_{0}=',num2str(d(end),'%6.5g')),'fontsize',30)
set(gca, 'FontSize', 15)
saveas(h,'eigenvec_0_paul.eps','epsc')
close all
h=figure(1);
semilogx(omega,eigenvec_anal(:,2),'b','linewidth',4)
hold on
semilogx(omega,-eigenvec_num(:,end-1),'r','linewidth',2)
lgd=legend('Analytical','Numerical');
set(lgd,'FontSize',20);
box off
grid off
xlim([0 20])
xlabel('\omega','fontsize',25)
title(strcat('\lambda_{1}=',num2str(d(end-1),'%6.5g')),'fontsize',30)
set(gca, 'FontSize', 15)
saveas(h,'eigenvec_1_paul.eps','epsc')
close all
h=figure(1);
semilogx(omega,eigenvec_anal(:,3),'b','linewidth',4)
hold on
semilogx(omega,eigenvec_num(:,end-2),'r','linewidth',2)
lgd=legend('Analytical','Numerical');
set(lgd,'FontSize',20);
box off
grid off
xlim([0 20])
xlabel('\omega','fontsize',25)
title(strcat('\lambda_{2}=',num2str(d(end-2),'%6.5g')),'fontsize',30)
set(gca, 'FontSize', 15)
saveas(h,'eigenvec_2_paul.eps','epsc')
close all
h=figure(1);
semilogx(omega,eigenvec_anal(:,4),'b','linewidth',4)
hold on
semilogx(omega,eigenvec_num(:,end-3),'r','linewidth',2)
lgd=legend('Analytical','Numerical');
set(lgd,'FontSize',20);
box off
grid off
xlim([0 20])
xlabel('\omega','fontsize',25)
title(strcat('\lambda_{3}=',num2str(d(end-3),'%6.5g')),'fontsize',30)
set(gca, 'FontSize', 15)
saveas(h,'eigenvec_3_paul.eps','epsc')
close all

% % Least squares error
% abs_error=zeros(4,1);
% for k=1:4
%     abs_error(k)=norm(eigenvec_anal(:,k)-eigenvec_num(:,end-k+1))*sqrt(domega);      
% end
% abs_error(2)=norm(eigenvec_anal(:,2)+eigenvec_num(:,end-2+1))*sqrt(domega);                
% fid=fopen('least_squares_error_eigenvectors_paul_5000.txt','w');   
% fprintf(fid,'%s\n','From the highest to the lowest eigenvalue:');
% for k=1:4
%     fprintf(fid,'%2.9g\n',abs_error(k));
% end
% fclose(fid);
% % Least squares error with more grid points (check sensitivity)
% % 10000 points
% omega=linspace(0,omega_num(end),10000);
% omega=omega';
% domega=omega(2)-omega(1);
% eigenvec_num=zeros(numel(omega),4);
% for k=1:4
%     eigenvec_num(:,end-k+1)=reconstruction_signal(V(:,end-k+1),domega_num,0,omega);
%     eigenvec_num(:,end-k+1)=eigenvec_num(:,end-k+1)./norm(eigenvec_num(:,end-k+1))./sqrt(domega);    
% end
% eigenvec_anal=paul_analytical_4_eigenvectors_frequency(omega);
% abs_error=zeros(4,1);
% for k=1:4
%     abs_error(k)=norm(eigenvec_anal(:,k)-eigenvec_num(:,end-k+1))*sqrt(domega);      
% end
% abs_error(2)=norm(eigenvec_anal(:,2)+eigenvec_num(:,end-2+1))*sqrt(domega);                
% fid=fopen('least_squares_error_eigenvectors_paul_10000.txt','w');   
% fprintf(fid,'%s\n','From the highest to the lowest eigenvalue:');
% for k=1:4
%     fprintf(fid,'%2.9g\n',abs_error(k));
% end
% fclose(fid);
% % 50000
% omega=linspace(0,omega_num(end),50000);
% omega=omega';
% domega=omega(2)-omega(1);
% eigenvec_num=zeros(numel(omega),4);
% for k=1:4
%     eigenvec_num(:,end-k+1)=reconstruction_signal(V(:,end-k+1),domega_num,0,omega);
%     eigenvec_num(:,end-k+1)=eigenvec_num(:,end-k+1)./norm(eigenvec_num(:,end-k+1))./sqrt(domega);    
% end
% eigenvec_anal=paul_analytical_4_eigenvectors_frequency(omega);
% abs_error=zeros(4,1);
% for k=1:4
%     abs_error(k)=norm(eigenvec_anal(:,k)-eigenvec_num(:,end-k+1))*sqrt(domega);      
% end
% abs_error(2)=norm(eigenvec_anal(:,2)+eigenvec_num(:,end-2+1))*sqrt(domega);                
% fid=fopen('least_squares_error_eigenvectors_paul_50000.txt','w');   
% fprintf(fid,'%s\n','From the highest to the lowest eigenvalue:');
% for k=1:4
%     fprintf(fid,'%2.9g\n',abs_error(k));
% end
% fclose(fid);



% Numerical values of the eigenvalues
myeig_anal=paul_analytical_eigenvalues(C,10);
disp([myeig_anal d(end:-1:end-9)])
fid=fopen('eigenvalues_paul_C_4.txt','w');    
fprintf(fid,'%s\n','Parameters:');
fprintf(fid,'%s %2.9g\n','eps_time',eps_time);
fprintf(fid,'%s %2.9g\n','eps_freq',eps_freq);
fprintf(fid,'%s\n','Analytical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',myeig_anal(k));
end
fprintf(fid,'%s\n','Numerical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',d(end-k+1));
end
fclose(fid);
C=30;
myeig_anal=paul_analytical_eigenvalues(C,10);
myeig_num=zeros(10,1);
for k=1:10
    myeig_num(k)=eigenvalue_frequency(V(:,end-k+1),C,omega_num,eps_time,eps_freq);
end
disp([myeig_anal myeig_num])
fid=fopen('eigenvalues_paul_C_30.txt','w');    
fprintf(fid,'%s\n','Parameters:');
fprintf(fid,'%s %2.9g\n','eps_time',eps_time);
fprintf(fid,'%s %2.9g\n','eps_freq',eps_freq);
fprintf(fid,'%s\n','Analytical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',myeig_anal(k));
end
fprintf(fid,'%s\n','Numerical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',myeig_num(k));
end
fclose(fid);
C=50;
myeig_anal=paul_analytical_eigenvalues(C,10);
myeig_num=zeros(10,1);
for k=1:10
    myeig_num(k)=eigenvalue_frequency(V(:,end-k+1),C,omega_num,eps_time,eps_freq);
end
disp([myeig_anal myeig_num])
fid=fopen('eigenvalues_paul_C_50.txt','w');    
fprintf(fid,'%s\n','Parameters:');
fprintf(fid,'%s %2.9g\n','eps_time',eps_time);
fprintf(fid,'%s %2.9g\n','eps_freq',eps_freq);
fprintf(fid,'%s\n','Analytical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',myeig_anal(k));
end
fprintf(fid,'%s\n','Numerical eigenvalues:');
for k=1:10
    fprintf(fid,'%2.9g\n',myeig_num(k));
end
fclose(fid);
