function Psi_t= Wavepacket_evolution(Psi_0, t, energy, eigenbasis)
% This function time-evolves the initial state Psi_0 using the eigenbasis
% matrix and the energy vector related to each eigenbasis vector.

Psi_t=zeros(length(energy),1); % Dimension of the wavepacket

for n=1:length(energy)   
    
    weight=conj(transpose(eigenbasis(:, n)))*Psi_0; % Projection on the eigenbasis  
    Psi_t=Psi_t+(weight*exp(-1i*t*energy(n)))*eigenbasis(:,n); 
end

