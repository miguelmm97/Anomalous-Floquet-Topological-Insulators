% ANALYSIS OF THE ANOMALOUS PHASE
% Here we see what is the nature of the edge states found inthe [-1,1]
% region of the phase diagram in the Kitagawa approach. First we are
% concerned with the distribution of the bare edge states, then we will
% explore the time dependence of wave-packets initialised on the edges of
% the system.

home
close all
clear all
addpath('functions')


%% Parameters
omega=4/3; % Frequency of the drive ([J])
delta=0; % Potential offset ([J])
lambda=3.4;% Degree of anisotropy
J=1; J_1=J; J_2=1; J_3=1; % Hoping amplitudes ([J])
termination='Zigzag';

% Time parameters
n_t=100;
T=2*pi/omega;
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
dt=t(2)-t(1);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;

% Geometric Parameters 
N_layers=26;
sublattice_layer=1:N_layers;

% Momentum range
kx_min=linspace(-pi,0,100);
kx_max=linspace(0,pi,100);
kx=[kx_min,kx_max];    


%% Band structure for the armchair geometry

for i=1:length(kx);   
   
    % Hamiltonian construction (we use the lattice site basis with the
    % order {cell, layer, sublattice}={(n,1,A),(n,2,A)...(n,m,A),(n,1,B)...(n,m,B)}
    
    [eigenstates, eigenvalues]=Floquet_rud_xperiodic(kx(i), N_layers, n_t, dt, J_1, J_2, J_3, delta, termination);
    bands(:,i)=diag(eigenvalues);

end
     

figure (2)
close all
box on
hold on
for j=1:2*N_layers
plot(kx, bands(j,:),'.b','MarkerSize',3)
end
%title(['\lambda=', num2str(lambda),'  \omega=',num2str(2*pi/T), ' \Delta /J =', num2str(delta)])
xlim([-pi pi])
ylim([-pi pi])
xlabel('$k_{\parallel}a$','FontSize',20,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',20,'FontWeight','bold','Interpret','latex')
ax = gca;
ax.FontSize =20; 
hold off


%% Wavefunctions of the edge states
% % We consider the edge states for a particular momentum k. We name the edge
% % states as edge_state_gap_position. As the spectrum is going to be
% % ordered, the pi states are always the states number 1 and end, and the 0
% % states are always the N and N+1. We use the lattice site basis.
% 
% k=1; % Units of 1/a
% 
% % Calculation of the  driven spectrum
% [eigenvectors, eigenvalues]=Time_evolution_edge_rud(k, N_layers, n_t, dt, J_1, J_2, J_3, delta); % Calculation of the floquet spectrum
% [A, B]=order_eigenvalues(eigenvalues, eigenvectors); % Order spectrum from lowest energy to highest
% bands=diag(A); states=B;
% 
% % Edge states
% % Because of the basis A,B we are using, the first N_layers entries of the
% % eigenvector will be localised in the A sublaticce, the others in the B
% % sublattice
% edge_state_pi_lower=states(:,1); 
% prob_density_pi_lower=abs(conj(edge_state_pi_lower).*edge_state_pi_lower);
% edge_state_pi_upper=states(:,end);
% prob_density_pi_upper=abs(conj(edge_state_pi_upper).*edge_state_pi_upper);
% edge_state_0_lower=states(:,N_layers);
% prob_density_0_lower=abs(conj(edge_state_0_lower).*edge_state_0_lower);
% edge_state_0_upper=states(:,N_layers+1);
% prob_density_0_upper=abs(conj(edge_state_0_upper).*edge_state_0_upper);
% 
% 
% figure(1)
% close all
% % subplot(2,1,1)
% % box on
% % hold on
% % plot(sublattice_layer(1:10), prob_density_pi_lower(1:10), '-.r', 'LineWidth', 3)
% % plot(sublattice_layer(1:10), prob_density_pi_lower(N_layers+1:N_layers+10), 'r', 'LineWidth', 3)
% % plot(sublattice_layer(end-10:end), prob_density_pi_upper(N_layers-10:N_layers), 'b', 'LineWidth', 3)
% % plot(sublattice_layer(end-10:end), prob_density_pi_upper(end-10:end), '-.b', 'LineWidth', 3)
% % text(13, 0.4, '$\pi$-gap','FontSize',20,'FontWeight','bold','Interpret','latex')
% % text(4, 0.3, 'B','FontSize',20,'FontWeight','bold','Interpret','latex')
% % text(28, 0.3, 'A','FontSize',20,'FontWeight','bold','Interpret','latex')
% % hold off
% % %title('$\pi$-gap','FontSize',20,'FontWeight','bold','Interpret','latex')
% % %legend( 'A (v>0)','B (v>0)','A (v<0)', 'B (v<0)','Location','Best')
% % xlabel('Layer Number','FontSize',20,'FontWeight','bold','Interpret','latex')
% % ylabel('$\vert \psi \vert^2$','FontSize',20,'FontWeight','bold','Interpret','latex')
% % ylim([0,0.5])
% % xlim([1 N_layers])
% % ax = gca;
% % ax.FontSize =20; 
% % 
% % subplot(2,1,2)
% box on
% hold on
% plot(sublattice_layer(end-10:end), prob_density_0_lower(N_layers-10:N_layers), 'b', 'LineWidth', 3)
% plot(sublattice_layer(end-10:end), prob_density_0_lower(end-10:end), '-.b', 'LineWidth', 3)
% plot(sublattice_layer(1:10), prob_density_0_upper(1:10), '-.r', 'LineWidth', 3)
% plot(sublattice_layer(1:10), prob_density_0_upper(N_layers+1:N_layers+10), 'r', 'LineWidth', 3)
% text(13, 0.4, '0-gap','FontSize',20,'FontWeight','bold','Interpret','latex')
% text(4, 0.3, 'B','FontSize',20,'FontWeight','bold','Interpret','latex')
% text(28, 0.3, 'A','FontSize',20,'FontWeight','bold','Interpret','latex')
% hold off
% %title('0-gap','FontSize',20,'FontWeight','bold','Interpret','latex')
% %legend( 'A (v>0)','B (v>0)','A (v<0)', 'B (v<0)','Location','Best')
% xlabel('Layer Number','FontSize',20,'FontWeight','bold','Interpret','latex')
% ylabel('$\vert \psi \vert^2$','FontSize',20,'FontWeight','bold','Interpret','latex')
% ylim([0,0.5])
% xlim([1 N_layers])
% ax = gca;
% ax.FontSize =20; 

% print -depsc edgechiralityloc[1,0].eps
