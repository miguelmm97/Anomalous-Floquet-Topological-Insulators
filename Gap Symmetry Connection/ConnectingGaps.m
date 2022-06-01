
%% SYMMETRY CONNECTION OF BOTH GAPS

% In this code we try to connect both gaps / branchcuts by applying certain
% symmetry operations to the Floquet Hamiltonian / Evolution operator

addpath('functions')

%% Parameters
clear all
home

% Parameters
omega=4.5; delta=2; lambda=3; T=2*pi/omega;
drive=1; % 0--> Rudner || 1--> Kitagawa
branchcut=1; % 0-> Branchcut at 0, 1-> Branchcut at pi quasienergy

J=1; J_1=J; J_2=J; J_3=J;

% Time parameters
n_t=100;
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
dt=t(2)-t(1);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;

% Momentum
kx=-3.1;ky=-0.56;


%% Connecting Branchcuts


switch drive
    case 0 % Rudner
        
        % 0-gap ( 0 branchcut)
        U=Time_evolution_bulk(delta, kx, ky, n_t, dt, J_1, J_2, J_3);
        HF_0=1i*logm(U); % Bare HF without the branchcut, default at 0
        
        % pi-gap ( 1 branchcut)
        [A, B]=eig(HF_0); [values, vectors]=order_eigenvalues(B, A);
        u_1=vectors(:,1);
        HF_pi= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)')... 
                               + real(values(2,2))* vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projectors)
               
    case 1 % Kitagawa
        
        % 0-gap ( 0 branchcut)
        U=Time_evolution_bulk_lambda(delta, kx, ky, n_t, dt, J_1, J_2, J_3, lambda);
        HF_0=1i*logm(U); % Bare HF without the branchcut, default at 0
        
        % pi-gap ( 1 branchcut)
        [A, B]=eig(HF_0); [values, vectors]=order_eigenvalues(B, A);
        u_1=vectors(:,1);
        HF_pi= real(values(1,1)+branchcut*2*pi)* vectors(:,1)*(vectors(:,1)')...
                               + real(values(2,2))* vectors(:,2)*(vectors(:,2)'); % Reconstructing HF (eigenvalue*projectors)
        
    otherwise
        disp('Error in the specified drive')
        
end