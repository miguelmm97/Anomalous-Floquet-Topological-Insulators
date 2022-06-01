function [eigenstates, eigenvalues]=Time_evolution_edge_rud(kx, N_layers, n_t, dt, J_1, J_2, J_3, delta)
% We calculate the time evolution of the armchair geometry at momentum kx,
% with N_layers layers, with period defined by n points in each driving step of
% length dt, and in the phase dictated by the three J_i couplings, the period,
% delta and lambda.

% Auxiliary vectors
main_diag=ones(N_layers,1);
sup_diag=ones(N_layers-1,1);
inf_diag=ones(N_layers-1,1);
diagonal=(delta/2)*ones(N_layers,1);

for s=1:N_layers
    if mod(s,2)==0
        main_diag(s)=exp(1i*kx);
    end
end   % Alternating Bloch Factor

% Off-diagonal matrices
gamma1=J_1*diag(main_diag,0);
gamma2=J_2*diag(inf_diag,-1);
gamma3=J_3*diag(sup_diag,1);

% Hamiltonians for each driving step
H1=[diag(diagonal),gamma1;
        conj(transpose(gamma1)), diag(-diagonal)];
    
    H2=[diag(diagonal),gamma2;
        conj(transpose(gamma2)),diag(-diagonal)];
    
    H3=[diag(diagonal),gamma3;
        conj(transpose(gamma3)),diag(-diagonal)];

% Evolution operators and Floquet Hamiltonian
U1=(expm(-1i*dt*H1))^n_t; U2=(expm(-1i*dt*H2))^n_t; U3=(expm(-1i*dt*H3))^n_t;
prod=U1*U2*U3; HF=1i*logm(prod);

[eigenstates, eigenvalues]=eig(HF);