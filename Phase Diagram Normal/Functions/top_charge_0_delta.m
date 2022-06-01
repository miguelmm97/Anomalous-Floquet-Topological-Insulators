function [q0, kx0, ky0]=top_charge_0_delta(cont_0, singularities_0, step, n_t, omega,J_1,J_2,J_3)


% Getting rid of redundant gap closings

control_0=compare(singularities_0(:,1),step);

% Calculation of winding numbers
cont=1;
if cont_0>1
for j=1:length(control_0(:,1))
    
    gap_aux=10;
    
    % No more groups of band touchings
    if isnan(control_0(j,1))==1;
        break
    end
    
    % Analysis of each group of touchings
    for p=1:length(control_0(1,:))
        
        if isnan(control_0(j,p))==1; % Group ends
            break
            
        else         
            ind=find( abs(singularities_0(:,1)-control_0(j,p))<0.0001); % Identification of gap
            
            if gap_aux>singularities_0(ind(1),5) % Condition to save gap, write ind(1) just in case there are 2 simultaneous gaps
                
                gap_aux=singularities_0(ind(1),5);
                indice=ind(1);                
            end
            
        end
    end
   
    % Topological charge of the group of band touchings
    delta=singularities_0(indice,1);
    kx=singularities_0(indice,2);
    ky=singularities_0(indice,3);   
    dt=singularities_0(indice,4);   
    h_0=0;
    q0(cont,:)=[ delta, S_matrix_delta(omega,kx,ky,delta,n_t,dt,h_0,J_1,J_2,J_3)];
    kx0(cont,:)=kx;
    ky0(cont,:)=ky;
    cont=cont+1;
end
end