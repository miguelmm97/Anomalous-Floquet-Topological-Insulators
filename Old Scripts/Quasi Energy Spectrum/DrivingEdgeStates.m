%% EDGE STATES

% We look at a strip geometry of 30-40 unitcells with armchair termination
% and momentum parallel to the zig-zag edge. In this case the armchair
% termination is for the x edge, and momentum parallel to the y edge is
% therefore k_x

%% 3-STEP PERIODIC DRIVING
clear all
home

Periods=[2*pi/3];

opts = odeset('RelTol',1e-2,'AbsTol',1e-2);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;


N_layers=32;
main_diag=ones(N_layers,1);
sup_diag=ones(N_layers-1,1);
inf_diag=ones(N_layers-1,1);
%% ARMCHAIR TERMINATION STATIC SYSTEM

kx_min=linspace(-pi,0,1000);
kx_max=linspace(0,pi,1000);
kx=[kx_min,kx_max];

for j=1:length(kx)
    
    % Construction of the hamiltonian
    H=zeros(2*N_layers,2*N_layers);
    gamma=zeros(N_layers,N_layers);
    
    k_par=kx(j);
    
for p=1:N_layers
        
        if mod(p,2)==0
            main_diag(p)=exp(1i*k_par);
        end
    end
    
    gamma=-diag(main_diag,0)-diag(sup_diag,1)-diag(inf_diag,-1);  
    H=[zeros(N_layers,N_layers),gamma;
        conj(gamma), zeros(N_layers,N_layers)];
    
    % Edge state dispersion
    E_edge(:,j)=eig(H);
end

figure (1)
box on
hold on
for j=1:2*N_layers
plot(kx,E_edge(j,:),'.b','MarkerSize',2)
end
hold off
title(['n=',num2str(N_layers)])
ax=gca;
ax = gca;
ax.FontSize = 16;
xlim([-pi pi])
ylim([-3 3])
xlabel('$k_{\parallel}a$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off

%% ARMCHAIR TERMINATION DRIVEN SYSTEM

kx_min=linspace(-pi,0,100);
kx_max=linspace(0,pi,100);
kx=[kx_min,kx_max];

Aux=diag(ones(2*N_layers,1));

main_diag=ones(N_layers,1);
sup_diag=ones(N_layers-1,1);
inf_diag=ones(N_layers-1,1);

g= @(t,y,A) -1i*A*y;
       
for j=1:length(Periods)
T=Periods(j);
% Driving step 1

for i=1:length(kx);   
   i 
    for s=1:N_layers
        if mod(s,2)==0
            main_diag(s)=exp(1i*kx(i));
        end
    end
    
    gamma=-diag(main_diag,0);
    
    H=[zeros(N_layers,N_layers),gamma;
        conj(gamma), zeros(N_layers,N_layers)];
    
    % Schrodinger's Equation for the 1st driving stage             
    for f=1:2*N_layers
        
    [t,y]=ode45(@(t,y) g(t,y,H) ,[0,T/3], Aux(:,f), opts);
    
    basis_1(:,f,i)=y(end,:);
    end

end



% Driving step 2
for i=1:length(kx);
    i
    gamma=-diag(inf_diag,-1);
    
    H=[zeros(N_layers,N_layers),gamma;
        conj(gamma), zeros(N_layers,N_layers)];
    
    % Schrodinger's Equation for the 1st driving stage    
    for f=1:2*N_layers
        
    [t,y2]=ode45(@(t,y2) g(t,y2,H) ,[T/3,2*T/3], basis_1(:,f,i), opts);
        
    basis_2(:,f,i)=y2(end,:);
    end

end

% Driving step 3
for i=1:length(kx);
    i
   gamma=-diag(sup_diag,1);
    
    H=[zeros(N_layers,N_layers),gamma;
        conj(gamma), zeros(N_layers,N_layers)];
    
    % Schrodinger's Equation for the 1st driving stage
    
    for f=1:2*N_layers
        
    [t,y3]=ode45(@(t,y3) g(t,y3,H) ,[2*T/3,T],basis_2(:,f,i),opts);
    
    basis_3(:,f,i)=y3(end,:);
    
    end
    
    for f=1:2*N_layers
     % Time evolution operator
     U=conj(Aux)*basis_3(:,:,i);
     U_diag=eig(U);
     % Quasi-Energy Bands
     eps_edge(f,i,j)=1i*log(U_diag(f))./Periods(j);
    end
end
     
end

%%

figure (2)
close all
box on
hold on
for j=1:2*N_layers
plot(kx,eps_edge(j,:,1).*Periods(1),'.b','MarkerSize',3)
end
title(['n=',num2str(N_layers),'  T=',num2str(Periods(1))])
ax=gca;
ax = gca;
ax.FontSize = 12;
xlim([-pi pi])
ylim([-pi pi])
xlabel('$k_{\parallel}a$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off


