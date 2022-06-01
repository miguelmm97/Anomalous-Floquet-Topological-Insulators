% TIGHT-BINDING GRAPHENE

clear all
home

% Reciprocal Lattice Parameters (Check Sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);


% Analytical dispersion
%E_an_Gamma_K=sqrt(1+(4*cos(3.*kx./2).*cos(sqrt(3).*ky./2))+4*(cos(sqrt(3).*ky./2)).^2);
%%

% M --> Gamma
kx=[a:-0.01:0];
ky=tan(30*pi/180)*kx;
for f=[1:length(kx)];
     
         H=[0, exp(-i.*ky(f))+(2.*exp(i.*ky(f)./2).*cos(sqrt(3).*kx(f)./2));
         exp(i.*ky(f))+(2.*exp(-i.*ky(f)./2).*cos(sqrt(3).*kx(f)./2)), 0];
              
    E_num_M_Gamma(f,:)=eig(H);
    E_an_M_Gamma(f)=sqrt(1+(4*cos(3.*ky(f)./2).*cos(sqrt(3).*kx(f)./2))+4*(cos(sqrt(3).*kx(f)./2)).^2);
    
end
% K --> M
kx=[l:-0.01:a];
ky=tan(120*pi/180)*kx-(l*tan(120*pi/180));
for f=[1:length(kx)];
   
       H=[0, exp(-i.*ky(f))+(2.*exp(i.*ky(f)./2).*cos(sqrt(3).*kx(f)./2));
         exp(i.*ky(f))+(2.*exp(-i.*ky(f)./2).*cos(sqrt(3).*kx(f)./2)), 0];
     
    % Numerical dispersion
    E_num_K_M(f,:)=eig(H);
    E_an_K_M(f)=sqrt(1+(4*cos(3.*ky(f)./2).*cos(sqrt(3).*kx(f)./2))+4*(cos(sqrt(3).*kx(f)./2)).^2);
    
end
% Gamma --> K
kx=[0:0.001:l];
ky=0;
for f=[1:length(kx)];
   
        H=[0, exp(-1i.*ky)+(2.*exp(1i.*ky./2).*cos(sqrt(3).*kx(f)./2));
         exp(1i.*ky)+(2.*exp(-1i.*ky./2).*cos(sqrt(3).*kx(f)./2)), 0];
     
    % Numerical dispersion
    E_num_Gamma_K(f,:)=eig(H);   
    E_an_Gamma_K(f)=sqrt(1+(4*cos(3.*ky./2).*cos(sqrt(3).*kx(f)./2))+4*(cos(sqrt(3).*kx(f)./2)).^2);
end


% Plotting parameters
K=l;
M=sqrt(b^2+(l-a)^2);
Gamma=2*pi/3;
k1=linspace(0,K,length(E_an_Gamma_K));
k2=linspace(K,K+M,length(E_an_K_M));
k3=linspace(K+M,K+M+Gamma,length(E_an_M_Gamma));

close all
figure(1)
box on
hold on
plot(k1,E_num_Gamma_K(:,2),'-r')
%plot(k1,E_an_Gamma_K,'--g')
plot(k2,E_num_K_M(:,2),'-r')
%plot(k2,E_an_K_M,'--g')
plot(k3,E_num_M_Gamma(:,2),'-r')
%plot(k3,E_an_M_Gamma,'--g')
%plot(zeros(length([0:3])),[0:3],'--k')
%legend('Numerical dispersion','Analytical dispersion')
ax=gca;
ax.XTick=[0,K,K+M,K+M+Gamma];
ax.XTickLabel={'\Gamma','K','M','\Gamma'}
ax = gca;
ax.FontSize = 16;
xlim([0 K+M+Gamma])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off
%print -djpeg -r600 bandstructure.png


%% Quasi Energy bands of the bulk, Other(?) High symmetry lines

% Gamma --> K2
kx=[0:0.01:s];
ky=tan(pi*60/180)*kx;
for j=1:length(kx);
    
     H=[0, exp(-i.*ky(j))+(2.*exp(i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2));
         exp(i.*ky(j))+(2.*exp(-i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2)), 0];
              
    eps_Gamma_K(j,:)=eig(H);
    
    
end
% K2 -->K
kx=[s:-0.01:-s];
ky=f*ones(1,length(kx));
for j=1:length(kx);
    H=[0, exp(-i.*ky(j))+(2.*exp(i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2));
         exp(i.*ky(j))+(2.*exp(-i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2)), 0];
              
    eps_K_K(j,:)=eig(H);
    
end
% K --> Gamma
kx=[-s:0.01:0];
ky=-tan(pi*60/180)*kx;
for j=1:length(kx);
    
    H=[0, exp(-i.*ky(j))+(2.*exp(i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2));
         exp(i.*ky(j))+(2.*exp(-i.*ky(j)./2).*cos(sqrt(3).*kx(j)./2)), 0];
              
    eps_K_Gamma(j,:)=eig(H);
    
    
end


% Plotting parameters
K=sqrt(s^2+f^2);
KK=2*s;
Gamma=K;
k1=linspace(0,K,length(eps_Gamma_K));
k2=linspace(K,K+KK,length(eps_K_K));
k3=linspace(K+KK,K+KK+Gamma,length(eps_K_Gamma));

figure(3)
box on
hold on
plot(k1,eps_Gamma_K,'.b')
plot(k2,eps_K_K,'.b')
plot(k3,eps_K_Gamma,'.b')
plot(K*ones(length([-3:3])),[-3:3],'--k')
plot(K+KK*ones(length([-3:3])),[-3:3],'--k')
hold off
ax=gca;
ax.XTick=[0,K,K+KK,K+KK+Gamma];
ax.XTickLabel={'\Gamma','K_1','K_2','\Gamma'}
ax = gca;
ax.FontSize = 16;
xlim([0 K+KK+Gamma])
ylim([-pi pi])
%title(['  \omega=',num2str(2*pi/T), ' \Delta /J =', num2str(delta)])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off


%% 3D BAND STRUCTURE
kx=[-4*l:0.01:4*l];
ky=[-4*l:0.01:4*l];

[X,Y]=meshgrid(kx);
E_an=sqrt(1+(4*cos(3.*X./2).*cos(sqrt(3).*Y./2))+4*(cos(sqrt(3).*Y./2)).^2);


for f=[1:length(kx)];
    for j=[1:length(ky)];
        
         H=[0, exp(-i.*kx(f))+(2.*exp(i.*kx(f)./2).*cos(sqrt(3).*ky(j)./2));
         exp(i.*kx(f))+(2.*exp(-i.*kx(f)./2).*cos(sqrt(3).*ky(j)./2)), 0];
       
       E_num(f,j,:)=eig(H);
        
    end
end
%%


kx=[-l:0.01:l];
ky=[-l:0.01:l];

J=1;
J1=J;
J2=0*J;
J3=0*J;
for f=[1:length(kx)];
    for j=[1:length(ky)];
        
         H=[0, -J1-J2*exp(-0.5*1i*(sqrt(3)*kx(f)+3*ky(j)))-J3*exp(-0.5*1i*(-sqrt(3)*kx(f)+3*ky(j)));
         -J1-J2*exp(0.5*1i*(sqrt(3)*kx(f)+3*ky(j)))-J3*exp(0.5*1i*(-sqrt(3)*kx(f)+3*ky(j))), 0];
       
       E_num(f,j,:)=eig(H);
        
    end
end

figure(2)
close all
hold on
box on
%mesh(E_an)
mesh(E_num(:,:,2))
mesh(E_num(:,:,1))

ax=gca;
ax.XTick=[0,242,484];
ax.XTickLabel={'-2.41','0','2.41','FontSize',18}
ax.YTick=[0,242,484];
ax.YTickLabel={'-2.41','0','2.41','FontSize',18}
ax = gca;
ax.FontSize = 16; 
% xlim([0 484])
% ylim([0 484])
xlabel('$k_x a$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$k_y a$','FontSize',18,'FontWeight','bold','Interpret','latex')
zlabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
%print -djpeg -r600 3Dbandstructure.png