%% EDGE STATES

% We look at a strip geometry of 30-40 unitcells with armchair termination
% and momentum parallel to the zig-zag edge. In this case the armchair
% termination is for the x edge, and momentum parallel to the y edge is
% therefore k_x

%% 3-STEP PERIODIC DRIVING
clear all
home

% Parameters
omega=0.;
T=2*pi/omega;
delta=2;
J_1=1;
J_2=1;
J_3=1.1;

% Time parameters
n_t=100;
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
dt=t(2)-t(1);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;


% Hamiltonin parameters
N_layers=32;
main_diag=ones(N_layers,1);
sup_diag=ones(N_layers-1,1);
inf_diag=ones(N_layers-1,1);
diagonal=(delta/2)*ones(N_layers,1);

%% ARMCHAIR TERMINATION DRIVEN SYSTEM

% Momentum range
kx_min=linspace(-pi,0,100);
kx_max=linspace(0,pi,100);
kx=[kx_min,kx_max];    

for i=1:length(kx);   
   
    % Hamiltonian construction
    for s=1:N_layers
        if mod(s,2)==0
            main_diag(s)=exp(1i*kx(i));
        end
    end    
    gamma1=J_1*diag(main_diag,0);
    gamma2=J_2*diag(inf_diag,-1);
    gamma3=J_3*diag(sup_diag,1);
    
    H1=[diag(diagonal),gamma1;
        conj(transpose(gamma1)), diag(-diagonal)];
   
    H2=[diag(diagonal),gamma2;
        conj(transpose(gamma2)),diag(-diagonal)];
    
    H3=[diag(diagonal),gamma3;
        conj(transpose(gamma3)),diag(-diagonal)];
  
       
     % Time evolution operator and quasi-energy
     U1=(expm(-1i*dt*H1))^n_t;
     U2=(expm(-1i*dt*H2))^n_t;
     U3=(expm(-1i*dt*H3))^n_t;
     prod=U1*U2*U3;
     U_diag=eig(prod);
     
     for f=1:2*N_layers    
     eps_edge(f,i)=1i*log(U_diag(f))./T;
    end
end
     

%%

figure (2)
close all
box on
hold on
for j=1:2*N_layers
plot(kx,eps_edge(j,:)*T,'.b','MarkerSize',3)
end
title(['n=',num2str(N_layers),'  \omega=',num2str(2*pi/T), ' \Delta /J =', num2str(delta)])
ax=gca;
ax = gca;
ax.FontSize = 12;
xlim([-pi pi])
ylim([-pi pi])
xlabel('$k_{\parallel}a$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off