% GREEN'S FUNCTIONS APPROACH TO GAP CHARACTERISATION

% In this code we compute the in-gap Green's functions for the Floquet topological
% phases following the thoery described by Robert.


home
close all
clear all
addpath('functions')

%% Parameters

% Parameters of the FTI
drive='Rudner';
omega=3; %(Units of J)
T=2*pi/omega;
delta=1; %(Units of J)
lambda=3.4;
J_1=1;
J_2=1;
J_3=1;
branchcut=0; % 0-> Branchcut at 0, 1-> Branchcut at pi quasienergy

% Time parameters
n_t=1000; % Resolution of the time vector used in the time evolution
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,3*T/3,n_t)]; % Time vector of the driv0
dt=t(2)-t(1); % Differential time step

% Definitions
kx=linspace(-pi, pi, 100); ky=0; dk=kx(2)-kx(1); % Perpendicular momentum along x
if branchcut==0
    eps=linspace(-pi, pi, 500); % Range of quasienergy BRANCHCUT 0 -PI--PI, BRANCHCUT 1 0--2PI
else
    eps=linspace(0, 2*pi, 500);
end



%% Matrix approach to GF
% We check that a topologically trivial hamiltonian with flatbands does
% indeed give no roots and the correct singularities in the GF by using
% direct matrix inverse.



for s=1:length(eps)   
    s
    % G(eps, kx, ky=0)
     for j=1:length(kx)    
         
         clear HF values vectors A B M
         
         if strcmp(drive, 'Rudner')==1
             
             U=Time_evolution_bulk(delta, kx(j),ky, n_t, dt, J_1, J_2, J_3);
             HF= [ zeros(2) ,      U;       % Doubled floquet hamiltonian
             ctranspose(U), zeros(2)];
            
             
         else
             U=Time_evolution_bulk_lambda(delta, kx(j),ky, n_t, dt, J_1, J_2, J_3, lambda);
             HF= [ zeros(2) ,      U;       % Doubled floquet hamiltonian
             ctranspose(U), zeros(2)];
         
         end
         
         M=(eps(s)*eye(4)-HF); % Denom in GF
         F(:,:,j)=inv(M)./(2*pi); % GF for eps_GF(s) at kx_GF (/2pi for later integration)
         
     end
     
     % In gap GF G(eps, ky=0)
     sum=0; % Initial value of the integration
     for t=1:length(kx)-1
     % Careful as is a matrix integration we cannot just use sum!!!!!!!!!!
     sum=sum+(F(:,:,t)+F(:,:,t+1))./2*dk; % Numerical integration of the matrix F for -pi<kx<pi
        
     end

    % GF eigenvalues
    eigenvalues(:,s)=eig(sum); % Eigenvalues of the Green's function
       
end



figure(2)
subplot(2,4,1:4)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(2,:),'.b')
plot(eps, eigenvalues(3,:),'.b')
plot(eps, eigenvalues(4,:),'.b')
plot(eps, 0*ones(1, length(eps)),'r') %, '    $\lambda$=', num2str(lambda)
title(['Greens function eigenvalues ','$\omega/J$=', num2str(omega), '    $\Delta/J=$', num2str(delta), '    $\lambda$=', num2str(lambda)],'FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
xlim([eps(1) eps(end)])
ax = gca;
ax.FontSize =15; 


subplot(2,4,5:6)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(2,:),'.b')
plot(eps, eigenvalues(3,:),'.b')
plot(eps, eigenvalues(4,:),'.b')
plot(eps, 0*ones(1, length(eps)),'r')
%title('0-gap','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
if branchcut==0 % Select zoom zone
    xlim([-0.4 0.4])
else
    xlim([2.5 4.5])
end
ax = gca;
ax.FontSize =15; 


subplot(2,4,7:8)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(2,:),'.b')
plot(eps, eigenvalues(3,:),'.b')
plot(eps, eigenvalues(4,:),'.b')
plot(eps, 0*ones(1, length(eps)),'r')
%title('$\pi$-gap','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
if branchcut==0 % Select zoom zone
    xlim([2.5 pi])
else
    xlim([-4.5 -2.5])
end
ax = gca;
ax.FontSize =15; 


set(gcf,'Position',[250 300 1000 600])

















