% Figures for chapter 3 of PhD thesis.

%%%%%%%%%%%%%%%%%%%%%%
%   Morlet Wavelet   %
%%%%%%%%%%%%%%%%%%%%%%

% Iso-contours: comparison ellipse vs exact solution - R=0.000001
w0=[2, 5.5, 10, 20];
R=0.000001;
first_guess=[0.00000000001 100000000000];
n_w0=numel(w0);
for k=1:n_w0
    [b_p_ap,b_p_am,b_m_ap,b_m_am,ap,am]=isocontour_morlet_exact(w0(k),R);
    [~,b_max,a_b_sq_max]=b_min_max(w0(k),R,first_guess);
    [~,~,a_center,a_spread]=a_min_max(w0(k),R,a_b_sq_max);       
    [x_el,y_el]=ellipse_a_b(b_max,a_center,a_spread);
    h=figure(k);
    hold on
    plot(b_p_ap,ap,'r','linewidth',3)
    plot(b_m_ap,ap,'b','linewidth',3)
    plot(b_p_am,am,'c','linewidth',3)
    plot(b_m_am,am,'g','linewidth',3)
    plot(x_el,y_el,'k--','linewidth',2)
    box off
    grid on
    xlim([-1.5*b_max 1.5*b_max])
    xlabel('b','fontsize',25)
    ylabel('a','fontsize',25)
    title(strcat('\omega_{0}=',num2str(w0(k))),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('iso_contours_morlet_R_small_w0_',num2str(w0(k)),'.eps'),'epsc')
    close all
end

% Iso-contours: comparison ellipse vs exact solution - R=0.5
w0=[0.5, 5.5, 10, 20];
R=0.5;
first_guess=[0.00001 100000];
n_w0=numel(w0);
for k=1:n_w0
    [b_p_ap,b_p_am,b_m_ap,b_m_am,ap,am]=isocontour_morlet_exact(w0(k),R);
    [~,b_max,a_b_sq_max]=b_min_max(w0(k),R,first_guess);
    [~,~,a_center,a_spread]=a_min_max(w0(k),R,a_b_sq_max);       
    [x_el,y_el]=ellipse_a_b(b_max,a_center,a_spread);
    h=figure(k);
    hold on
    plot(b_p_ap,ap,'r','linewidth',3)
    plot(b_m_ap,ap,'b','linewidth',3)
    plot(b_p_am,am,'c','linewidth',3)
    plot(b_m_am,am,'g','linewidth',3)
    plot(x_el,y_el,'k--','linewidth',2)
    box off
    grid on
    xlim([-1.5*b_max 1.5*b_max])
    xlabel('b','fontsize',25)
    ylabel('a','fontsize',25)
    title(strcat('\omega_{0}=',num2str(w0(k))),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('iso_contours_morlet_R_medium_w0_',num2str(w0(k)),'.eps'),'epsc')
    close all
end

% Iso-contours: comparison ellipse vs exact solution - R=0.8
w0=[0.0000001, 5.5, 10, 20];
R=0.8;
first_guess=[0.00001 100000];
n_w0=numel(w0);
for k=1:n_w0
    [b_p_ap,b_p_am,b_m_ap,b_m_am,ap,am]=isocontour_morlet_exact(w0(k),R);
    [~,b_max,a_b_sq_max]=b_min_max(w0(k),R,first_guess);
    [~,~,a_center,a_spread]=a_min_max(w0(k),R,a_b_sq_max);       
    [x_el,y_el]=ellipse_a_b(b_max,a_center,a_spread);
    h=figure(k);
    hold on
    plot(b_p_ap,ap,'r','linewidth',3)
    plot(b_m_ap,ap,'b','linewidth',3)
    plot(b_p_am,am,'c','linewidth',3)
    plot(b_m_am,am,'g','linewidth',3)
    plot(x_el,y_el,'k--','linewidth',2)
    box off
    grid on
    xlim([-1.5*b_max 1.5*b_max])
    xlabel('b','fontsize',25)
    ylabel('a','fontsize',25)
    title(strcat('\omega_{0}=',num2str(w0(k))),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('iso_contours_morlet_R_big_w0_',num2str(w0(k)),'.eps'),'epsc')
    close all
end

% Iso-contours for various values of R. w0 is fixed
w0=6;
R=[0.001 0.01 0.1 0.5 0.9];
first_guess=2;
n_R=numel(R);
h=figure(1);
hold on
for k=1:n_R
    [b_min,b_max,a_b_sq_max]=b_min_max(w0,R(k),first_guess);
    [a_min,a_max,a_center,a_spread]=a_min_max(w0,R(k),a_b_sq_max);       
    [x_el,y_el]=ellipse_a_b(b_max,a_center,a_spread);
    plot(x_el,y_el,'k','linewidth',1)
    plot(b_min,a_b_sq_max,'r.','markersize',20)
    plot(b_max,a_b_sq_max,'r.','markersize',20)
    text(b_max*0.1,a_max+0.05,strcat('R=',num2str(R(k))),'fontsize',15)
end
box off
grid on
axis([-6.5 6.5 0 3])
xlabel('b','fontsize',20)
ylabel('a','fontsize',20)
set(gca, 'FontSize', 15)
saveas(h,'iso_contours_morlet_w0_fixed.eps','epsc')
close all

% c_psi
myw0=5.5:0.001:50;
mycpsi=zeros(numel(myw0),1);
for k=1:numel(myw0)
    mycpsi(k)=c_phi(myw0(k));
end
h=figure(1);
plot(myw0,mycpsi,'k','linewidth',2)
box off
grid off
xlim([0 50])
xlabel('\omega_0','fontsize',20)
ylabel('c_{\psi}','fontsize',20)
set(gca, 'FontSize', 15)
saveas(h,'c_psi.eps','epsc')
close all

% First 9 eigenfunctions - w0=6 and R=0.01 and 0.001
R=0.01;
w0=6;
eps=0.000001;   
n_eigenvec=5000; % number of frequency grid points for all the eigenvectors
eigenvec_num=zeros(n_eigenvec,9);
[omega_num1,d1,V1]=integrated_projection_simple_integral_exact_contour_frequency(w0,R,eps,1);
close all
domega_num1=omega_num1(2)-omega_num1(1);
omega1=linspace(0,omega_num1(end),n_eigenvec);
omega1=omega1';
domega1=omega1(2)-omega1(1);
mysign=ones(9,1);
mysign(3)=-1;
mysign(4)=-1;
mysign(8)=-1;
for k=1:9
    eigenvec_num(:,end-k+1)=reconstruction_signal(V1(:,end-k+1),domega_num1,omega_num1(1),omega1);
    eigenvec_num(:,end-k+1)=eigenvec_num(:,end-k+1)./norm(eigenvec_num(:,end-k+1))./sqrt(domega1);
    h=figure(1);
    hold on
    plot(omega1,mysign(k).*eigenvec_num(:,end-k+1),'k','linewidth',2)
    box off
    grid off
    xlim([0 20])
    ylim([-0.8 0.8])
    xlabel('\omega','fontsize',25)
    title(strcat('\lambda_{',num2str(k-1),'}=',num2str(d1(end-k+1),5)),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('eigenvec_',num2str(k-1),'_morlet_w0_6_R_0,01.eps'),'epsc')
    close all
end 

% First 3 eigenfunctions - w0=15 and R=0.001
R=0.001;
w0=15;
eps=0.000001;   
n_eigenvec=5000; % number of frequency grid points for all the eigenvectors
eigenvec_num=zeros(n_eigenvec,3);
[omega_num1,d1,V1]=integrated_projection_simple_integral_exact_contour_frequency(w0,R,eps,1);
close all
domega_num1=omega_num1(2)-omega_num1(1);
omega1=linspace(0,omega_num1(end),n_eigenvec);
omega1=omega1';
domega1=omega1(2)-omega1(1);
mysign=ones(3,1);
for k=1:3
    eigenvec_num(:,end-k+1)=reconstruction_signal(V1(:,end-k+1),domega_num1,omega_num1(1),omega1);
    eigenvec_num(:,end-k+1)=eigenvec_num(:,end-k+1)./norm(eigenvec_num(:,end-k+1))./sqrt(domega1);
    h=figure(1);
    hold on
    plot(omega1,mysign(k).*eigenvec_num(:,end-k+1),'k','linewidth',2)
    box off
    grid off
    xlim([10 30])
    ylim([-0.8 0.8])
    xlabel('\omega','fontsize',25)
    title(strcat('\lambda_{',num2str(k-1),'}=',num2str(d1(end-k+1),5)),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('eigenvec_',num2str(k-1),'_morlet_w0_15_R_0,001.eps'),'epsc')
    close all
end

% Asymptotic invariance properties and Hermite functions - work with the
% 4th eigenvector
R=[0.1 0.001];
colorplot=['b','k'];
w0=[6 15 30 100];
eps=0.000001;   
n_eigenvec=10000;
eigenvec_num=zeros(n_eigenvec,numel(R));
eigenval=zeros(numel(R),1);
omega=zeros(n_eigenvec,numel(R));
for l=1:numel(w0)
    for k=1:numel(R)
        [omega_num,d,V]=integrated_projection_simple_integral_exact_contour_frequency(w0(l),R(k),eps,1);
        close all
        domega_num=omega_num(2)-omega_num(1);
        omega(:,k)=linspace(omega_num(1),omega_num(end),n_eigenvec);
        domega=omega(2,k)-omega(1,k);
        eigenvec_num(:,k)=reconstruction_signal(V(:,end-3),domega_num,omega_num(1),omega(:,k));
        eigenvec_num(:,k)=eigenvec_num(:,k)./norm(eigenvec_num(:,k))./sqrt(domega); 
        eigenval(k)=d(end-3);
    end
    h=figure(1);
    hold on
    if l==2
        plot(omega(:,1),eigenvec_num(:,1),colorplot(1),'linewidth',1)
        plot(omega(:,2),eigenvec_num(:,2),colorplot(2),'linewidth',1)
    else
        for k=1:numel(R)
            plot(omega(:,k),(-1)^(k-1).*eigenvec_num(:,k),colorplot(k),'linewidth',1)
        end
    end
    plot(omega(:,2),-hermite_function(3,omega(:,2)-w0(l)),'r','linewidth',1)
    box off
    grid off
    xlim([w0(l)-6 w0(l)+6])
    ylim([-0.8 0.8])
    xlabel('\omega','fontsize',25)
    lgd=legend(strcat('R=',num2str(R(1)),' - \lambda_3=',num2str(eigenval(1))),strcat('R=',num2str(R(2)),' - \lambda_3=',num2str(eigenval(2))),'Hermite fun. \phi_3(\omega-\omega_0)');
    set(lgd,'fontsize',15)
    title(strcat('\omega_0=',num2str(w0(l))),'fontsize',30)
    set(gca, 'FontSize', 15)
    saveas(h,strcat('morlet_approx_invariance_4th_eigenvec_properties_w0_',num2str(w0(l)),'.eps'),'epsc')
    close all
end 

% % Asymptotic invariance properties and Hermite functions - work with the
% % 8th eigenvector
% R=[0.1 0.001];
% colorplot=['b','k'];
% w0=[6 15 30 100];
% eps=0.000001;   
% n_eigenvec=10000;
% eigenvec_num=zeros(n_eigenvec,numel(R));
% eigenval=zeros(numel(R),1);
% omega=zeros(n_eigenvec,numel(R));
% mysign=ones(numel(R),numel(w0));
% mysign(2,2)=-1;
% mysign(2,3)=-1;
% mysign(2,4)=-1;
% for l=1:numel(w0)
%     for k=1:numel(R)
%         [omega_num,d,V]=integrated_projection_simple_integral_exact_contour_frequency(w0(l),R(k),eps,1);
%         close all
%         domega_num=omega_num(2)-omega_num(1);
%         omega(:,k)=linspace(omega_num(1),omega_num(end),n_eigenvec);
%         domega=omega(2,k)-omega(1,k);
%         eigenvec_num(:,k)=reconstruction_signal(V(:,end-7),domega_num,omega_num(1),omega(:,k));
%         eigenvec_num(:,k)=eigenvec_num(:,k)./norm(eigenvec_num(:,k))./sqrt(domega); 
%         eigenval(k)=d(end-7);
%     end
%     h=figure(1);
%     hold on
%     for k=1:numel(R)
%         plot(omega(:,k),mysign(k,l).*eigenvec_num(:,k),colorplot(k),'linewidth',1)
%     end
%     plot(omega(:,2),-hermite_function(7,omega(:,2)-w0(l)),'r','linewidth',1)
%     box off
%     grid off
%     xlim([max(0,w0(l)-8) w0(l)+8])
%     xlabel('\omega','fontsize',20)
%     legend(strcat('R=',num2str(R(1)),' - \lambda_7=',num2str(eigenval(1))),strcat('R=',num2str(R(2)),' - \lambda_7=',num2str(eigenval(2))),'Hermite fun. \phi_7(\omega-\omega_0)','fontsize',20)
%     title(strcat('\omega_0=',num2str(w0(l))),'fontsize',20)
%     set(gca, 'FontSize', 15)
%     saveas(h,strcat('morlet_approx_invariance_properties_8th_eigenvec_w0_',num2str(w0(l)),'.eps'),'epsc')
%     close all
% end
