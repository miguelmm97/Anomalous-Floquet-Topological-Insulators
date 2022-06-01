%% GREEN'S FUNCTIONS APPROACH TO GAP CHARACTERISATION

% In this code we compute the in-gap Green's functions for the Floquet topological
% phases following the theory described by Robert.


home
close all
clear all
addpath('functions')

% Parameters

% Parameters of the FTI
drive='Rudner';
omega=4/3; %(Units of J)
T=2*pi/omega;
delta=0.5; %(Units of J)
lambda=3.4;
J_1=1; J_2=1; J_3=1;
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
%HF=zeros(2,2,length(kx)); % Floquet hamiltonian definition
h0=zeros(1,length(kx)); % Floquet vector definition
hx=zeros(1,length(kx)); % Floquet vector definition
hy=zeros(1,length(kx)); % Floquet vector definition
hz=zeros(1,length(kx)); % Floquet vector definition

sigmax=[0, 1; 1, 0]; % Pauli matrix x
sigmay=[0, -i; i, 0]; % Pauli matrix y
sigmaz=[1, 0; 0, -1]; % Pauli matrix z

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
             
             U=Time_evolution_bulk_2T(delta, kx(j),ky, n_t, dt, J_1, J_2, J_3);
             HF=1i*logm(U); % Bare HF without the branchcut, default at 0
             % Taking the specified branchcut
             [A, B]=eig(HF); [values, vectors]=order_eigenvalues(B, A);
             u_1(:,j)=vectors(:,1);
             HF= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)') + real(values(2,2))* vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projectors)
             
         else
             U=Time_evolution_bulk_lambda_2T(delta, kx(j),ky, n_t, dt, J_1, J_2, J_3, lambda);
             HF=1i*logm(U);
             % Taking the specified branchcut
             [A, B]=eig(HF); [values, vectors]=order_eigenvalues(B, A);
             HF= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)') + real(values(2,2))*vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projecto
             u_1(:,j)=vectors(:,1);
         end
         
         M=(eps(s)*eye(2)-HF); % Denom in GF
         F(:,:,j)=inv(M)./(2*pi); % GF for eps_GF(s) at kx_GF (/2pi for later integration)
         
     end
     
     % In gap GF G(eps, ky=0)
     sum=0; % Initial value of the integration
     for t=1:length(kx)-1
     % Carefull as is a matrix integration we cannot just use sum!!!!!!!!!!
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
plot(eps, 0*ones(1, length(eps)),'r')
%title('0-gap','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
if branchcut==0 % Select zoom zone
    xlim([-1.3 1.3])
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

%% Pauli approach to GF

% Floquet Hamiltonian
for i=1:length(kx)
    
    clear HF values vectors A B 
    if strcmp(drive, 'Rudner')==1
        
        U=Time_evolution_bulk_2T(delta, kx(i),ky, n_t, dt, J_1, J_2, J_3);
        HF=1i*logm(U); % Bare HF without the branchcut, default at 0
        % Taking the specified branchcut
        [A, B]=eig(HF); [values, vectors]=order_eigenvalues(B, A);
        u_1(:,i)=vectors(:,1);
        HF= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)') + real(values(2,2))* vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projectors)
        
    else
        U=Time_evolution_bulk_lambda_2T(delta, kx(i),ky, n_t, dt, J_1, J_2, J_3, lambda);
        HF=1i*logm(U);
        % Taking the specified branchcut
        [A, B]=eig(HF); [values, vectors]=order_eigenvalues(B, A);
        HF= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)') + real(values(2,2))*vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projecto
        u_1(:,i)=vectors(:,1);
    end
    
    % Pauli decomposition
    hx(i)=real(HF(2,1)); h0(i)=(HF(1,1)+HF(2,2))/2;
    hy(i)=imag(HF(2,1)); hz(i)=(HF(1,1)-HF(2,2))/2;
    h2(i)=hx(i)^2+hy(i)^2+hz(i)^2;
end

%Calculation of the Green's function
for s=1:length(eps)
    s
    % Integrands for g, gx, gy, gz
    for j=1:length(kx)
    fx(j)=hx(j)/(((eps(s)-h0(j))^2-h2(j))*2*pi);
    fy(j)=hy(j)/(((eps(s)-h0(j))^2-h2(j))*2*pi);
    fz(j)=hz(j)/(((eps(s)-h0(j))^2-h2(j))*2*pi);
    f(j)=(eps(s)-h0(j))/(((eps(s)-h0(j))^2-h2(j))*2*pi);
    end
    
    % Integrals for g, g0, gx, gy, gz
    Gx(s)=sum(fx(1:length(kx)-1)+ fx(2:length(kx)))/2*dk;
    Gy(s)=sum(fy(1:length(kx)-1)+ fy(2:length(kx)))/2*dk;
    Gz(s)=sum(fz(1:length(kx)-1)+ fz(2:length(kx)))/2*dk;
    G(s)=sum(f(1:length(kx)-1)+ f(2:length(kx)))/2*dk;
    
    % Complete Green's function
    G_tot=G(s)*eye(2)+Gx(s)*sigmax+Gy(s)*sigmay+Gz(s)*sigmaz;
    eigenvalues(:,s)=eig(G_tot); % Eigenvalues of the Green's function
    
end



figure(2)
subplot(2,4,1:4)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(2,:),'.b')
plot(eps, 0*ones(1, length(eps)),'r') %, '    $\lambda$=', num2str(lambda)
title(['Greens function eigenvalues ','$\omega/J$=', num2str(omega), '    $\Delta/J=$', num2str(delta), '    $\lambda$=', num2str(lambda)],'FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET/J$','FontSize',15,'FontWeight','bold','Interpret','latex')
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
plot(eps, 0*ones(1, length(eps)),'r')
%title('0-gap','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET/J$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
if branchcut==0 % Select zoom zone
    xlim([-2 2])
else
    xlim([2 4.5])
end
ax = gca;
ax.FontSize =15; 


subplot(2,4,7:8)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(2,:),'.b')
plot(eps, 0*ones(1, length(eps)),'r')
%title('$\pi$-gap','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET/J$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
if branchcut==0 % Select zoom zone
    xlim([2.5 pi])
else
    xlim([5.9 2*pi])
end
ax = gca;
ax.FontSize =15; 

set(gcf,'Position',[250 300 1000 600])
