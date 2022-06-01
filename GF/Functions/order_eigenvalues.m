function [C,D]=order_eigenvalues(A,B)
% This functions uses the eigenvalue matrix A and orders its eigenvalues
% (in the main diagonal) from lowest to greatest, and afterwards orders the
% eigenvector matrix B in the corresponding order.

cont=2; % Loop counter

aux_1=real(diag(A)); final=length(aux_1); % Need to impose reality condition bc if not the minimum is not defined
eigenvalues=zeros(1, final); eigenvectors=zeros(final, final);

% First iteration
[eigenvalues(1), idx]=min(aux_1); % Minimum eigenvalue goes first
aux_1(idx)=NaN; % Delete used entry
eigenvectors(:,1)=B(:, idx); % Minimum eigenvector goes first

while cont<=final    
   clear idx
   [eigenvalues(cont), idx]=min(aux_1);
   aux_1(idx)=NaN;   
   eigenvectors(:, cont)=B(:, idx);
   cont=cont+1;     
end

C=diag(eigenvalues); % Matrix form
D=eigenvectors;
