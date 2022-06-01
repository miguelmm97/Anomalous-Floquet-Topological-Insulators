%% PHASE DIAGRAM FLOQUET GRAPHENE (omega)

clear all
home
addpath('Functions')

% IMPORTANT COMMENTS: When looking at the anisotropic case, it is necessary
% to look at all the high symmetry lines and set a gap tolerance of around
% 0.1 at least, when looking at the isotropic case it is enough to look at
% the high symmetry points, and set a tolerance of around 0.01

%% Parameters
% Drive Family Physical Parameters
step=0.01; % Resolution of the frequency vector
omega=4; %  Frequency (units of J)
T=2*pi./omega; % Period 
delta=[0:0.1:2]; % Potential offset vector
lambda=[2:step:5]; % Degree of anisotropy (lambda=1 identified with high frequency static regime)
J=1; J_1=J; J_2=J; J_3=J; % Every quantity in terms of J


% Computational parameters
tol=0.01; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vehomctor used in the time evolution
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
dt=t(2)-t(1); % Time resolution of the drive

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% High Symmetry Lines 
kx_gamma_K=0;   % From the Gamma point to K
ky_gamma_K=0;
kx_KK=l;        % From K to K'
ky_KK=0;
kx_K_gamma=s;   % From K to the Gamma point
ky_K_gamma=f;
% 
% kx_gamma_K=[s:-0.1:0];   % From K2 to the Gamma point
% ky_gamma_K=tan(60*pi/180)*kx_gamma_K;
% kx_K_gamma=[0:0.1:l];    % From the Gamma point to K1
% ky_K_gamma=0*kx_K_gamma;
% kx_KK=[l:-0.1:s];      % From K1 to K2
% ky_KK=tan(120*pi/180)*kx_KK-(l*tan(120*pi/180));
%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.

for z=1:length(delta)
    
    disp(delta(z))
    W_0=0;  % Initial 0 winding number at omega(1)
    W_pi=1; % Initial pi winding number at omega(2)    
    
    % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    cont_pi=1; % Number of times a singularity pi is detected
    clear singularities_0 singularities_pi q0 qpi idx kx_aux ky_aux kx_idx ky_idx
       
    for j=1:length(lambda)
       
        % Time evolution
        for i=1:length(kx_gamma_K);
            
            U=Time_evolution_bulk_lambda(delta(z),kx_gamma_K(i),ky_gamma_K(i),n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_Gamma_K_1(i)=HF_diag(1)/J; % Quasi energy bands
            eps_Gamma_K_2(i)=HF_diag(2)/J; % Quasi energy bands
            
        end
        for i=1:length(kx_KK);
            
            U=Time_evolution_bulk_lambda(delta(z), kx_KK(i), ky_KK(i), n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_K_1(i)= HF_diag(1)/J; % Quasi energy bands
            eps_K_K_2(i)= HF_diag(2)/J; % Quasi energy bands
            
            
        end
        for i=1:length(kx_K_gamma);
            
            U=Time_evolution_bulk_lambda(delta(z), kx_K_gamma(i), ky_K_gamma(i), n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_Gamma_1(i)= HF_diag(1)/J; % Quasi energy bands
            eps_K_Gamma_2(i)= HF_diag(2)/J; % Quasi energy bands
            
        end
        
        % Gap calculation
        gap_0=min([abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]);
        gap_pi=min(pi*ones(1,length([abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]))-[abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]);               
        if  gap_0<tol
            idx = find([abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]==gap_0); % Position of the gap
            kx_aux=[kx_gamma_K, kx_KK, kx_K_gamma];
            ky_aux=[ky_gamma_K, ky_KK, ky_K_gamma];
            kx_idx=kx_aux(idx(1)); % x-momentum of the gap, idx(1) in case there are 2 simultaneous gaps
            ky_idx=ky_aux(idx(1)); % y-momentum of the gap
            
            % Point in 3D parameter space where the band touching occurs
            singularities_0(cont_0,:) = [lambda(j), kx_idx, ky_idx, dt, gap_0]; % We save dt for the S matrix
            cont_0=cont_0+1; % 1 extra singularity has been detected
        end       
        if  gap_pi<tol
            
            idx = find(pi*ones(1,length([abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]))-[abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]==gap_pi);
            kx_aux=[kx_gamma_K, kx_KK, kx_K_gamma];
            ky_aux=[ky_gamma_K, ky_KK, ky_K_gamma];
            kx_idx=kx_aux(idx(1)); % x-momentum of the gap
            ky_idx=ky_aux(idx(1)); % y-momentum of the gap
            
            E_aux_1=[abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]; % Calculation of h_0 in order to decompose on Pauli matrices
            E_aux_2=[abs(eps_Gamma_K_2),abs(eps_K_K_2),abs(eps_K_Gamma_2)]; % Calculation of h_0 in order to decompose on Pauli matrices
            h_0=-(E_aux_1(idx(1))+E_aux_2(idx(1)))/2; % Calculation of h_0 in order to decompose on Pauli matrices
            
            
            % Point in 3D parameter space where the band touching occurs
            singularities_pi(cont_pi,:) = [lambda(j), kx_idx, ky_idx, dt, h_0]; % We save dt for the S matrix
            cont_pi=cont_pi+1; % 1 extra singularity has been detected
        end
        
    end
    
    % Singularities and charges    
    if cont_0>1
        q0=top_charge_0_lambda(cont_0, singularities_0, step, n_t, delta(z),J_1,J_2,J_3); % Assign a charge to each singularity
        
    else
        q0=0;
    end   
    if cont_pi>1
        qpi=top_charge_pi_lambda(cont_pi, singularities_pi, step, n_t, delta(z),J_1,J_2,J_3); % Assign a charge to each singularity
      
    else
        qpi=0;
    end  

    % Phase diagram
    for j=1:length(lambda)
        
        for i=1:length(q0(:,1))
            if lambda(j)==q0(i,1)
                W_0=W_0-q0(i,2); % Adding the topological charge at each singularity
            end
        end
        for i=1:length(qpi(:,1))
            if lambda(j)==qpi(i,1)
                W_pi=W_pi-qpi(i,2);  % Adding the topological charge at each singularity
            end
        end
        
        phase_diag(z,j,:)=[delta(z), lambda(j), W_0, W_pi];
        
    end
    
end

%% Figures

figure(2)
close all
hold on
box on

for z=1:length(delta)
    
    for j=1:1:length(lambda)
        
        if phase_diag(z,j,3)==0 && phase_diag(z,j,4)==0 ;
            plot( lambda(j), delta(z), '.b','MarkerSize',5)
            
        else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==0 ;
                plot(lambda(j), delta(z), '.g','MarkerSize',5)
                
            else if phase_diag(z,j,3)==0 && phase_diag(z,j,4)==1 ;
                    plot(lambda(j), delta(z), '.','Color',[0 0.8 0.2],'MarkerSize',5)
                    
                    
                else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==1 || phase_diag(z,j,3)==-1 && phase_diag(z,j,4)==-1 ;
                        plot(lambda(j), delta(z), '.r','MarkerSize',5)
                        
                    else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==-1 ;
                            plot(lambda(j), delta(z), '.k','MarkerSize',5)
                            
                            else if  phase_diag(z,j,3)==-1 && phase_diag(z,j,4)==1;
                            plot(lambda(j), delta(z), '.k','MarkerSize',5)
                            
                        else if phase_diag(z,j,3)==-1 && phase_diag(z,j,4)==0 ;
                                plot(lambda(j), delta(z), '.y','MarkerSize',5)
                                
                                
                            else if  phase_diag(z,j,3)==0 && phase_diag(z,j,4)==-1;
                                    plot(lambda(j), delta(z), '.c','MarkerSize',5)
                                    
                                else
                                    plot(lambda(j), delta(z), '.m','MarkerSize',5)
                                end
                            end
                        end
                    end
                end
                end
            end
        end
    end
end

title(['Phase Diagram for w=', num2str(omega)])
xlabel('$\lambda$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')














