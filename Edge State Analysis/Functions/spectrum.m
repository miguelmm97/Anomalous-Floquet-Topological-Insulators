function [energy, eigenstates] = spectrum(H)
% Diagonalisation of H

[states, values]=eig(H);                                                % Diagonal Hamiltonian
[eigenvalues, eigenstates]=order_eigenvalues(values, states);           % Ordered spectrum
energy=diag(eigenvalues);                                               % Energies