function [HF, H1, H2, H3] = Floquet_kit_strip(delta, lambda, J1, J2, J3, Nx, Ny, position_matx, T, t0)
% This function constructs the Floquet hamiltonian of a finite strip provided the
% couplings, delta and lambda in the Kitagawa drive.

% Step wise hamiltonians
H1=Finite_Hamiltonian( delta, lambda*J1, J2, J3, Nx, Ny, position_matx); % 1st step Kitagawa
H2=Finite_Hamiltonian( delta, J1, lambda*J2, J3, Nx, Ny, position_matx); % 2nd step Kitagawa
H3=Finite_Hamiltonian( delta, J1, J2, lambda*J3, Nx, Ny, position_matx); % 3rd step Kitagawa

% Stroboscopic time evolution (Floquet)
if t0 <= T/3
    t1 = (T/3) - t0; t2 = T/3; t3 = T/3; t4 = t0;
    U1 = expm(-1i*t1*H1);
    U2 = expm(-1i*t2*H2);
    U3 = expm(-1i*t3*H3);
    U4 = expm(-1i*t4*H1);
else if T/3 < t0 <= 2*T/3
    t1 = (2*T/3) - t0; t2 = T/3; t3 = T/3; t4 = t0 - T/3;
    U1 = expm(-1i*t1*H2);
    U2 = expm(-1i*t2*H3);
    U3 = expm(-1i*t3*H1);
    U4 = expm(-1i*t4*H2);
    else
    t1 = T - t0; t2 = T/3; t3 = T/3; t4 = t0 - (2*T/3);
    U1 = expm(-1i*t1*H3);
    U2 = expm(-1i*t2*H1);
    U3 = expm(-1i*t3*H2);
    U4 = expm(-1i*t4*H3);
    end
    
end
    
U=U4*U3*U2*U1; 
HF=(1i/T)*logm(U); % Floquet Hamiltonian  
