function position_matx = position_matrix(Nx, Ny)
% This funciton constructs the position matrix of a finite strip provided the
% the number of layers, where layers are counted vertically 
% moving to the right of the strip. Cells are refer to as the hexagons of
% the honeycomb grid.

Ns=Nx*Ny/2;     % Number of states
y_sep=sqrt(3)/2; % Separation between y layers in units of a=1
position_matx=zeros(Ns,3); % (n_state, xpos, ypos)

% For every site (state) ...
for n=1:Ns
    
    position_matx(n,1)=n; % State number   
    a=ceil(n/Ny);         % Selection of the layer (sweep over sites in vertical direction)
    cell=ceil(a/2)-1;     % Number of the cell the state is in
    
    
     % Two left x layers of the cell
    if mod(a,2)~=0
               
        % Inner left layer
        if a*Ny-n<=(Ny-1)/2 
            if a~=1
                position_matx(n,2)=1/2+3*cell; % Normal x position
            else
                position_matx(n,2)=1/2; % First cell x position
            end
            position_matx(n,3)=2*(a*Ny-n)*y_sep; % y position
           
        % Outer left layer
        else 
            if a~=1
                position_matx(n,2)=3*cell; % Normal x position
            else
                position_matx(n,2)=0;      % First cell x position
            end
            position_matx(n,3)=(2*((Ny-1)/2-Ny+a*Ny-n)+1)*y_sep; % y position
        end
        
        
    % Two right layers of the cell    
    else  
        
         % Outer right layer
        if a*Ny-n<(Ny-1)/2
            if a~=2
                position_matx(n,2)=2+3*cell; % Normal x position
            else
                position_matx(n,2)=2; % First cell x position
            end
            position_matx(n,3)=(2*(a*Ny-n)+1)*y_sep; % y position
            
        % Inner right layer    
        else 
            if a~=2
                position_matx(n,2)=3/2+3*cell; % Normal x position
            else
                position_matx(n,2)=3/2; % First cell x position
            end
            position_matx(n,3)=(2*((Ny-1)/2-Ny+a*Ny-n)+2)*y_sep; % y position
        end
        
    end
      
end 