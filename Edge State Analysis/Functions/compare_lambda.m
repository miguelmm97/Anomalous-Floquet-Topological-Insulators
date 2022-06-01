function [vec]=compare_lambda(frec,steps)
% This function compares the elements of the frequency vector, and groups
% them together if they are contiguous. This is measured by the step
% argument, two elements are contiguous if they differ by the step size.
% The function groups them in the rows of a frec x frec matrix, returning
% NaN values in the rest of the elements


% Variables
lambda=frec;
step=steps;
vec=ones(length(lambda),length(lambda)).*NaN; 

% Parameters
tol=0.00001;
cont1=1;
cont2=1;

for j=1:length(lambda)    
    
    if j==1
        vec(cont1,cont2)=lambda(j);
        cont2=cont2+1;
        
      else if abs(lambda(j)-control)<tol % Very important, never compare two decimals
            vec(cont1,cont2)=lambda(j);
            cont2=cont2+1;
         else
            cont1=cont1+1;
            cont2=1;
            vec(cont1,cont2)=lambda(j);
            cont2=cont2+1;            
      end
    end
    
    control=lambda(j)+step; % Contiguous value
end


