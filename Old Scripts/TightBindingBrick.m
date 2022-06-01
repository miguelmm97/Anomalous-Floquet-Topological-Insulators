%% DIAGONALISATION TIGHT BINDING BRICK LATTICE A,B SUBLATTICES

clear all
home

kx=[0:0.01:pi]; % in units of lattice parameter
ky=[0:0.01:pi]; % in units of lattice parameter

%% Gamma --> M

ky=0;
% Analytical dispersion
E_an_Gamma_M=sqrt(1+(4*cos(kx).*cos(ky))+4*(cos(kx)).^2);

% Numerical dispersion
E_num_Gamma_M=zeros(length(kx),2);

for f=[1:length(kx)];
     
        H=[0, exp(-i*ky)+2*cos(kx(f));
           exp(i*ky)+2*cos(kx(f)), 0];
              
    E_num_Gamma_M(f,:)=eig(H);
    
end

%% M --> X

kx=[pi:-0.01:pi/2];
ky=kx-pi;

E_num_M_X=zeros(length(kx),2);
E_an_M_X=zeros(length(kx),1);

% Numerical dispersion
for f=[1:length(kx)];
   
        H=[0, exp(-i*ky(f))+2*cos(kx(f));
           exp(i*ky(f))+2*cos(kx(f)), 0];
     
    % Numerical dispersion
    E_num_M_X(f,:)=eig(H);
    % Analytical dispersion
    E_an_M_X(f)=sqrt(1+(4*cos(kx(f)).*cos(ky(f)))+4*(cos(kx(f))).^2);
    
end

%% X--> Gamma

kx=[pi/2:-0.01:0];
ky=-kx;

for f=[1:length(kx)];
   
        H=[0, exp(-i*ky(f))+2*cos(kx(f));
           exp(i*ky(f))+2*cos(kx(f)), 0];
     
    % Numerical dispersion
    E_num_X_Gamma(f,:)=eig(H);
    % Analytical dispersion
    E_an_X_Gamma(f)=sqrt(1+(4*cos(kx(f)).*cos(ky(f)))+4*(cos(kx(f))).^2);
    
end

%% BAND STRUCTURE
kx=[0:0.01:pi];
%k2=linspace(pi,pi+pi/sqrt(2),length(E_an_M_X));
k3=linspace(-pi/sqrt(2),0,length(E_an_X_Gamma));

close all
figure(1)
box on
hold on
plot(kx,E_num_Gamma_M(:,2),'-r')
plot(kx,E_an_Gamma_M,'--g')
%plot(k2,E_an_M_X','--g')
%plot(k2,E_num_M_X(:,2),'--b')
plot(k3,E_num_X_Gamma(:,2),'-r')
plot(k3,E_an_X_Gamma,'--g')
plot(zeros(length([0:3])),[0:3],'--k')

legend('Numerical dispersion','Analytical dispersion')

ax=gca;
ax.XTick=[-pi/sqrt(2),0,pi];
ax.XTickLabel={'X','\Gamma','M'}
ax = gca;
ax.FontSize = 16;
xlim([-pi/sqrt(2) pi])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off
print -djpeg -r600 bandstructure.png


% 3D BAND STRUCTURE
kx=[-pi:0.01:pi];
ky=[-pi:0.01:pi];

[X,Y]=meshgrid(kx);
E_an=sqrt(1+(4*cos(X).*cos(Y))+4*(cos(X)).^2);

for f=[1:length(kx)];
    for j=[1:length(ky)];
        
         H=[0, exp(-i*ky(j))+2*cos(kx(f));
           exp(i*ky(j))+2*cos(kx(f)), 0];
       
       E_num(f,j,:)=eig(H);
        
    end
end


figure(2)
hold on
box on
mesh(E_num(:,:,2))
mesh(E_num(:,:,1))

ax=gca;
ax.XTick=[0,314,628];
ax.XTickLabel={'-\pi','0','\pi','FontSize',18}
ax.YTick=[0,314,628];
ax.YTickLabel={'-\pi','0','\pi','FontSize',18}
ax = gca;
ax.FontSize = 16; 
xlim([0 628])
ylim([0 628])
xlabel('$k_x a$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$k_y a$','FontSize',18,'FontWeight','bold','Interpret','latex')
zlabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
print -djpeg -r600 3Dbandstructure.png


