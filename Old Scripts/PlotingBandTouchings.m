%% PLOTING THE BAND TOUCHINGS

clear all
home

% Drive Family Parameters
omega=[0.9]; % units of J
periods=2*pi./omega;
delta=0.2;
tol=0.01;
W_0=0;
W_pi=0;

% Time parameters
n_t=1000;
change=0.5;

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% High Symmetry Lines 
kx_gamma_K=[0:0.01:l];
ky_gamma_K=0*kx_gamma_K;
kx_KK=[l:-0.01:s];
ky_KK=tan(120*pi/180)*kx_KK-(l*tan(120*pi/180));
kx_K_gamma=[s:-0.01:0];
ky_K_gamma=tan(60*pi/180)*kx_K_gamma;

% Hamiltonians

H1=@(delta,kx,ky)   [delta/2, exp(1i*ky);
                    exp(-1i*ky), -delta/2];    
                
H2=@(delta,kx,ky)  [delta/2, exp(-1i*0.5*(sqrt(3)*kx+ky));
                   exp(1i*0.5*(sqrt(3)*kx+ky)), -delta/2];    
               
H3=@(delta,kx,ky)  [delta/2, exp(-1i*0.5*(-sqrt(3)*kx+ky));
                    exp(1i*0.5*(-sqrt(3)*kx+ky)), -delta/2];

% S Matrix
cont_0=1;
cont_pi=1;

%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.
for j=1:length(omega)
j
% Time parameters
T=periods(j);   
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
dt=t(2)-t(1);

% Time evolution
for i=1:length(kx_gamma_K);
    
    % Time evolution operator
     U1=(expm(-1i*dt*H1(delta,kx_gamma_K(i),ky_gamma_K(i))))^n_t;
     U2=(expm(-1i*dt*H2(delta,kx_gamma_K(i),ky_gamma_K(i))))^n_t;
     U3=(expm(-1i*dt*H3(delta,kx_gamma_K(i),ky_gamma_K(i))))^n_t;
     prod=U1*U2*U3;
     
     HF=1i*logm(prod);%T
     HF_diag=eig(HF);
     
     % Quasi-Energy Bands
     eps_Gamma_K_1(i)=HF_diag(1);
     eps_Gamma_K_2(i)=HF_diag(2);
end
for i=1:length(kx_KK);
      
     % Time Evolution operator
     U1=(expm(-1i*dt*H1(delta,kx_KK(i),ky_KK(i))))^n_t;
     U2=(expm(-1i*dt*H2(delta,kx_KK(i),ky_KK(i))))^n_t;
     U3=(expm(-1i*dt*H3(delta,kx_KK(i),ky_KK(i))))^n_t;
     prod=U1*U2*U3;
     
     HF=1i*logm(prod);
     HF_diag=eig(HF);
     
     % Quasi-Energy Bands
     eps_K_K_1(i)=HF_diag(1);
     eps_K_K_2(i)=HF_diag(2);
     
end
for i=1:length(kx_K_gamma);
    
     % Time evolution operator
     U1=(expm(-1i*dt*H1(delta,kx_K_gamma(i),ky_K_gamma(i))))^n_t;
     U2=(expm(-1i*dt*H2(delta,kx_K_gamma(i),ky_K_gamma(i))))^n_t;
     U3=(expm(-1i*dt*H3(delta,kx_K_gamma(i),ky_K_gamma(i))))^n_t;
     prod=U1*U2*U3;
     
     HF=1i*logm(prod);
     HF_diag=eig(HF);
     
     % Quasi-Energy Bands
     eps_K_Gamma_1(i)=HF_diag(1);
     eps_K_Gamma_2(i)=HF_diag(2);
     
end


end


%% QUASI-ENERGY BAND STRUCTURE

K=l;
KK=sqrt(f^2+(l-s)^2);
Gamma=2*pi/3;

k1=linspace(0,K,length(eps_Gamma_K_1));
k2=linspace(K,K+KK,length(eps_K_K_1));
k3=linspace(K+KK,K+KK+Gamma,length(eps_K_Gamma_1));

T=periods;

close all
figure(1)
box on
hold on
plot(k1,eps_Gamma_K_1,'b.',k1,eps_Gamma_K_2,'b.')
plot(k2,eps_K_K_1,'b.',k2,eps_K_K_2,'b.')
plot(k3,eps_K_Gamma_1,'b.',k3,eps_K_Gamma_2,'b.')
plot(K*ones(length([-3:3])),[-3:3],'--k')
plot(K+KK*ones(length([-3:3])),[-3:3],'--k')
hold off

%legend('Numerical dispersion')

ax=gca;
ax.XTick=[0,K,K+KK,K+KK+Gamma];
ax.XTickLabel={'\Gamma','K_{+}','K_{-}','\Gamma'}
ax = gca;
ax.FontSize = 16;
xlim([0 K+KK+Gamma])
ylim([-pi pi])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off