function [eigenstates, eigenvalues]=Floquet_kit_xperiodic(kx, N_layers, n_t, dt, J_1, J_2, J_3, delta, lambda, termination)
% We calculate the time evolution of the armchair geometry at momentum kx,
% with N_layers layers, with period defined by n points in each driving step of
% length dt, and in the phase dictated by the three J_i couplings, the period,
% delta and lambda.

switch termination
    
    case 'Armchair'
       
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
        H1=[diag(diagonal),(lambda*gamma1)+gamma2+gamma3;
            conj(transpose((lambda*gamma1)+gamma2+gamma3)), diag(-diagonal)];
        H2=[diag(diagonal),gamma1+(lambda*gamma2)+gamma3;
            conj(transpose(gamma1+(lambda*gamma2)+gamma3)),diag(-diagonal)];
        H3=[diag(diagonal),gamma1+gamma2+(lambda*gamma3);
            conj(transpose(gamma1+gamma2+(lambda*gamma3))),diag(-diagonal)];
        
        % Evolution operators and Floquet Hamiltonian
         U1=expm(-1i*H1*T/3); U1=expm(-1i*H2*T/3); U3=U1=expm(-1i*H3*T/3);
         prod=U3*U2*U1; HF=1i*logm(prod);  
        [eigenstates, eigenvalues]=spectrum(HF);
        
        
    case 'Zigzag' % (Symmetrix zigzag edge)
        % Auxiliary vectors
        diagonal1=ones(floor(N_layers/2),1);
        diagonal2=ones(floor(N_layers/2),1);
        diagonal3=ones(floor(N_layers/2),1);
        sup_diag=ones(floor(N_layers/2)-1,1);              
            
        for s=1:floor(N_layers/2)
            if mod(s,2)==0
                diagonal2(s)=exp(-1i*kx*sqrt(3));
            else
                diagonal3(s)=exp(1i*kx*sqrt(3));
            end
        end   % Alternating Bloch Factors       
        
        % Off-diagonal matrices
        gamma1 = diag(sup_diag, 1).*J1;
        gamma2 = diag(diagona2,0).*J2;
        gamma3 = diag(diagona3,0).*J3;
              
        % Hamiltonians for each driving step
        H1=[0.5*delta*diag(diagonal), lambda*gamma1 + gamma2 + gamma3 ;
            conj(transpose(lambda*gamma1 + gamma2 + gamma3)), -0.5*delta*diag(diagonal)];
        H1=[0.5*delta*diag(diagonal), gamma1 + lambda*gamma2 + gamma3 ;
            conj(transpose(gamma1 + lambda*gamma2 + gamma3)), -0.5*delta*diag(diagonal)];
        H1=[0.5*delta*diag(diagonal), gamma1 + gamma2 + lambda*gamma3 ;
            conj(transpose(gamma1 + gamma2 + lambda*gamma3)), -0.5*delta*diag(diagonal)];

        
        % Evolution operators and Floquet Hamiltonian
        U1=expm(-1i*H1*T/3); U1=expm(-1i*H2*T/3); U3=U1=expm(-1i*H3*T/3);
        prod=U3*U2*U1; HF=1i*logm(prod);  
        [eigenstates, eigenvalues]=spectrum(HF);
        
        
        
end



