function H=Finite_Hamiltonian(delta, J1, J2, J3, Nx, Ny, position_matx)
% This funciton constructs the hamiltonian of a finite strip provided the
% couplings, the number of layers, where layers are counted vertically 
% moving to the right of the strip starting at 0, and the position of each state. 
% Remember position matrix =(n_state, xpos, ypos) and states are labelled
% from up to down in the same layer.

Ns=Nx*Ny/2;      % Total number of states
H=zeros(Ns,Ns);  % Hamiltonian
y_sep=sqrt(3)/2; % Separation between y layers in units of a=1

% For every site(state)
for i=1:Ns
    
    %%%%%%%%% Sublattice offset
    if mod(position_matx(i,2), 3)==0 || mod(position_matx(i,2)-3/2, 3)==0 
        H(i, i)=delta/2;  % A sublattice, located at the 1st and 3rd layers of the hexagons
    else
        H(i, i)=-delta/2; % B sublattice, located at the 2nd and 4th layers of the hexagons
    end
      
    
    %%%%%%%%% Couplings
    b=ceil(i/Ny);     % Selection of each group of Ny elements
    cell=ceil(b/2)-1; % Cell selection
       
    %%%% Edges in y direction 
    % Layer x=0 
    if i<=(Ny-1)/2
        H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny-1)/2>
        H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
        
    % Layer x=end
    else if i>Ns-(Ny-1)/2
            H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
            H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny-1)/2>   
            
        %%%% Inner layers sweeping in x direction
        % Odd groups of layers
        else if mod(b,2)~=0 
                
                %%%% Edges in x direction                
                % Layer y=0
                if position_matx(i,3)==0
                    H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                    H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                    
                % Layer y=end
                else if position_matx(i,3)==y_sep*(Ny-1)
                        H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
                        H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                        
                    %%%% Bulk outer layer 
                    else if position_matx(i,2)==3*cell; 
                            H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
                            H(i, i-(Ny-1)/2)=J1;  % Coupling of state |i> with |i-(Ny-1)/2>
                            
                        %%%%% Bulk inner layer
                        else
                            H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                            H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                            
                        end
                    end
                end
                
                
            % Even groups of layers
            else if mod(b,2)==0 
                    
                    %%%% Edges in x direction  
                    % Layer y=0
                    if position_matx(i,3)==0
                        H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny+1)/2>
                        H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i-(Ny+1)/2>
                        
                    % Layer y=end    
                    else if position_matx(i,3)==y_sep*(Ny-1)
                            H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i-(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
                            
                        %%%% Bulk long column   
                        else if position_matx(i,2)==3/2+3*cell
                                H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                                H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                                H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny+1)/2>
                                
                            %%%% Bulk short column    
                            else
                                H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
                                H(i, i+(Ny-1)/2)=J1;  % Coupling of state |i> with |i+(Ny-1)/2>
                                H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny-1)/2>
                                
                            end
                        end
                    end
                end
            end
        end
    end
end 

