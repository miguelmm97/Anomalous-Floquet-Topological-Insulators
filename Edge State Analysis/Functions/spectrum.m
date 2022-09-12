function [energy, eigenstates] = spectrum(H)
% Diagonalisation of H

[states, values]=eig(H);  % Diagonal Hamiltonian
% In case the eigenstates do not form a complete basis we use schur algorythm
if any(abs(1-diag(states'*states))>1e-6) || max(max(abs(states'*states-eye(length(H)))))>1e-6
    [states, t_aux] = schur(H);
    if max(abs(imag(t_aux)))>1e-7
        error('imaginary t_aux')
    end
    values = real(t_aux);
end

[eigenvalues, eigenstates]=order_eigenvalues(values, states);           % Ordered spectrum
energy=diag(eigenvalues);                                               % Energies
