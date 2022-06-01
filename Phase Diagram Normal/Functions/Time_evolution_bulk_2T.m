function  [prod]=Time_evolution_bulk_2T(delta,kx,ky,n_t,dt,J_1,J_2,J_3)
% For a given potentia offset, delta, point in the BZ, kx ky, and with a
% resolution in time given by n_t and dt, and the three hoping amplitudes (in units of J) 
% we calculate the Floquet Hamiltonian of the bulk (If we want to calculate for another system or
% for the edge system, just replace the Hamiltonians with the proper ones.

% Hamiltonians
H1=@(delta,kx,ky)   [delta/2, J_1*exp(1i*ky);  % Hamiltonian with only first hoping
                    J_1*exp(-1i*ky), -delta/2];    
                
H2=@(delta,kx,ky)  [delta/2, J_2*exp(-1i*0.5*(sqrt(3)*kx+ky)); % Hamiltonian with only second hoping
                   J_2*exp(1i*0.5*(sqrt(3)*kx+ky)), -delta/2];    
               
H3=@(delta,kx,ky)  [delta/2, J_3*exp(-1i*0.5*(-sqrt(3)*kx+ky)); % Hamiltonian with only third hoping
                    J_3*exp(1i*0.5*(-sqrt(3)*kx+ky)), -delta/2];
               


 % Time Evolution using the evolution operator
 U1=(expm(-1i*dt*H1(delta,kx,ky)))^n_t;
 U2=(expm(-1i*dt*H2(delta,kx,ky)))^n_t;
 U3=(expm(-1i*dt*H3(delta,kx,ky)))^n_t;
 prod=U1*U2*U3*U1*U2*U3;