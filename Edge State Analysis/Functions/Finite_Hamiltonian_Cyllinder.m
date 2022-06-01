function H=Finite_Hamiltonian_Cyllinder(delta, J1, J2, J3, Nx, Ny, position_matx)
% This funciton constructs the hamiltonian of a finite strip provided the
% couplings, the number of layers and the position of each state

Ns=Nx*Ny/2; % Total number of states
H=zeros(Ns,Ns); % Hamiltonian
y_sep=sqrt(3)/2; % Separation between y layers so that every site is equidistant from the others in units of a=1

for i=1:Ns
    
    % Sublattice offset, REMEMBER WE START COUNTING 0 LAYERS
    if mod( position_matx(i,2), 3)==0 || mod( position_matx(i,2)-3/2, 3)==0
        H(i, i)=delta/2;
    else
        H(i, i)=-delta/2;
    end
      
    % Coupling entries, REMEMBER WE START COUNTING 0 LAYERS
    b=ceil(i/Ny); % Selection of each group of Ny elements
    cell=ceil(b/2)-1;
       
    % Layer x=0 
    if i<=(Ny-1)/2
        H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny-1)/2>
        H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
        H(i, Ns-(Ny-1)/2+i)=J1; % Coupling of the first layer with the last layer
        H(Ns-(Ny-1)/2+i, i)=J1; % Coupling of the last layer with the first layer
    % Layer x=end
    else if i>Ns-(Ny-1)/2
            H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
            H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny-1)/2>
           
            
        else if mod(b,2)~=0 % Odd groups
                
                % Layer y=0
                if position_matx(i,3)==0
                    H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                    H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                % Layer y=end
                else if position_matx(i,3)==y_sep*(Ny-1)
                        H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
                        H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                    % Bulk outer layer 
                    else if position_matx(i,2)==3*cell; %position_matx(i,2)==2*b-2
                            H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
                            H(i, i-(Ny-1)/2)=J1;  % Coupling of state |i> with |i-(Ny-1)/2>
                        % Bulk inner layer
                        else
                            H(i, i-(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                            H(i, i-(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                            
                        end
                    end
                end
                
            else if mod(b,2)==0 % Even groups
                    
                    % Layer y=0
                    if position_matx(i,3)==0
                        H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i+(Ny+1)/2>
                        H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i-(Ny+1)/2>
                    % Layer y=end    
                    else if position_matx(i,3)==y_sep*(Ny-1)
                            H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i-(Ny-1)/2>
                            H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i+(Ny+1)/2>
                        % Bulk long column   
                        else if position_matx(i,2)==3/2+3*cell
                                H(i, i-(Ny+1)/2)=J1;  % Coupling of state |i> with |i+(Ny+1)/2>
                                H(i, i+(Ny+1)/2)=J2;  % Coupling of state |i> with |i-(Ny+1)/2>
                                H(i, i+(Ny-1)/2)=J3;  % Coupling of state |i> with |i-(Ny+1)/2>
                            % Bulk short column    
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
end % Bare Hamiltonian Construction

