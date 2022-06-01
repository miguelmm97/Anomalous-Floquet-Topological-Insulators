function [HF, H1, H2, H3] = Floquet_kit_cyl(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t)
% This function constructs the Floquet hamiltonian of a finite strip provided the
% couplings, delta and lambda in the Kitagawa drive.

% Step wise hamiltonians
H1=Finite_Hamiltonian_Cyllinder( delta, lambda*J_1, J_2, J_3, Nx, Ny, position_matx); % 1st step Kitagawa
H2=Finite_Hamiltonian_Cyllinder( delta, J_1, lambda*J_2, J_3, Nx, Ny, position_matx); % 2nd step Kitagawa
H3=Finite_Hamiltonian_Cyllinder( delta, J_1, J_2, lambda*J_3, Nx, Ny, position_matx); % 3rd step Kitagawa

% Stroboscopic time evolution (Floquet)
U1=(expm(-1i*dt*H1))^n_t; % Evolution operator 3
U2=(expm(-1i*dt*H2))^n_t; % Evolution operator 3
U3=(expm(-1i*dt*H3))^n_t; % Evolution operator 3
U=U3*U2*U1; 
HF=1i*logm(U); % Floquet Hamiltonian  