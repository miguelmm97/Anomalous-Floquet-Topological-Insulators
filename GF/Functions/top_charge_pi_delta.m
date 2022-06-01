function  [qpi, kxpi, kypi]=top_charge_pi_delta( cont_pi, singularities_pi, step, n_t, omega,J_1,J_2,J_3)


% Getting rid of redundant gap closings

control_pi=compare(singularities_pi(:,1),step);

% Calculation of winding numbers

cont=1;
if cont_pi>1
for j=1:length(control_pi(:,1))
      
      gap_aux=10;
      
       
    % No more groups of band touchings
    if isnan(control_pi(j,1))==1;
        break
    end
       
    
    % Analysis of each group of touchings
    for p=1:length(control_pi(1,:))
        
        if isnan(control_pi(j,p))==1;   % Group ends        
            break
        else
            
        ind=find( abs(singularities_pi(:,1)-control_pi(j,p))<0.0001); % Identification of gap
        
        if gap_aux>abs(pi-singularities_pi(ind(1),5)) % Condition to save gap, write ind(1) just in case there are 2 simultaneous gaps
            gap_aux=abs(pi-singularities_pi(ind(1),5));
            indice=ind(1);
        end
        end
    end
    
    % Topological charge of the group of band touchings
    delta=singularities_pi(indice,1);
    kx=singularities_pi(indice,2);
    ky=singularities_pi(indice,3);   
    dt=singularities_pi(indice,4);   
    h_0=singularities_pi(indice,5);
    
    qpi(cont,:)=[ delta, S_matrix_delta(omega,kx,ky,delta,n_t,dt,h_0,J_1,J_2,J_3)];
    kxpi(cont,:)=kx;
    kypi(cont,:)=ky;
    cont=cont+1;
end
end 