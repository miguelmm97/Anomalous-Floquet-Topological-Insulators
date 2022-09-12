function evolution_matx = Edge_evolution(prob_density, position_matx, Ny, delta_x, edge)
% This function gives back a matrix containing the evolution of the edge
% probability in real space. Given the full evolution in prob_density(n, T)
% and the position matrix, we get evolution_matx(T, x) where x is in
% batches of delta_x

N = length(prob_density(:,1)); Nt = length(prob_density(1,:));              % Number of states, number of evolution periods
y_sep=sqrt(3)/2; L = max(position_matx(:,2));                               % Separtion between y lyers and length of the strip in x dir
Nx = floor(L / delta_x);                                                    % Number of discretised x points
evolution_matx = zeros(Nt, Nx);                                             % Declare evolution matrix

for n=1:N 
    switch edge
        
        case 'Upper'
            if position_matx(n,3)>y_sep*(Ny-2)-0.1
                indx = floor(position_matx(n, 2)/delta_x) + 1;
                evolution_matx(:, indx)= prob_density(n,:);
            end
            
        case 'Lower'
            if position_matx(n,3)<y_sep + 0.1
               indx = floor(position_matx(n, 2)/delta_x) + 1;
               evolution_matx(:, indx)= prob_density(n,:);
            end
            
        otherwise
            disp('Unknown edge!')
    end
end


% cont=1;
% for n=1:Ns_real 
%     switch edge
%         case 'Upper'
%             if position_matx(n,3)>y_sep*(Ny-2)-0.1
%                 evolution_matx(:,cont)=prob_density(n,:);
%                 cont=cont+1;
%             end
%         case 'Lower'
%             if position_matx(n,3)<y_sep + 0.1
%                 evolution_matx(:,cont)=prob_density(n,:);
%                 cont=cont+1;
%             end
%             
%         otherwise
%             disp('Unknown edge!')
%     end
% end