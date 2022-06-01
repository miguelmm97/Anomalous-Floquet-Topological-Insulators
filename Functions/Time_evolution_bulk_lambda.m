function  [prod]=Time_evolution_bulk_lambda(delta,kx,ky,n_t,dt,J_1,J_2,J_3, lambda)
% For a given potentia offset, delta, point in the BZ, kx ky, and with a
% resolution in time given by n_t and dt, and the three hoping amplitudes (in units of J) 
% we calculate the Floquet Hamiltonian of the bulk (If we want to calculate for another system or
% for the edge system, just replace the Hamiltonians with the proper ones.)

% Hamiltonians


H1=@(delta,kx,ky)   [delta/2, (lambda*J_1)*exp(1i*ky)+J_2*exp(-1i*0.5*(sqrt(3)*kx+ky))+J_3*exp(-1i*0.5*(-sqrt(3)*kx+ky));  % Hamiltonian first step
                    (lambda*J_1)*exp(-1i*ky)+J_2*exp(1i*0.5*(sqrt(3)*kx+ky))+J_3*ex p(1i*0.5*(-sqrt(3)*kx+ky)), -delta/2];    
                
H2=@(delta,kx,ky)  [delta/2, J_1*exp(1i*ky)+(lambda*J_2)*exp(-1i*0.5*(sqrt(3)*kx+ky))+J_3*exp(-1i*0.5*(-sqrt(3)*kx+ky));  % Hamiltonian second step
                    J_1*exp(-1i*ky)+(lambda*J_2)*exp(1i*0.5*(sqrt(3)*kx+ky))+J_3*exp(1i*0.5*(-sqrt(3)*kx+ky)), -delta/2];    
               
H3=@(delta,kx,ky)  [delta/2, J_1*exp(1i*ky)+J_2*exp(-1i*0.5*(sqrt(3)*kx+ky))+(lambda*J_3)*exp(-1i*0.5*(-sqrt(3)*kx+ky));  % Hamiltonian third step
                    J_1*exp(-1i*ky)+J_2*exp(1i*0.5*(sqrt(3)*kx+ky))+(lambda*J_3)*exp(1i*0.5*(-sqrt(3)*kx+ky)), -delta/2];    
               


 % Time Evolution using the evolution operator
 U1=(expm(-1i*dt*H1(delta,kx,ky)))^n_t;
 U2=(expm(-1i*dt*H2(delta,kx,ky)))^n_t;
 U3=(expm(-1i*dt*H3(delta,kx,ky)))^n_t;
 prod=U3*U2*U1;