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
step=0.001; % Resolution of the frequency vector
omega=[7:-0.2:4.5]; %  Frequency (units of J) % Period 
delta=[2]; % Potential offset vector
lambda=[2:step:5]; % Degree of anisotropy (lambda=1 identified with high frequency static regime)
J=1; J_1=J; J_2=J; J_3=J; % Every quantity in terms of J


% Computational parameters
tol=0.01; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vehomctor used in the time evolution

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

for z=1:length(omega)
    
    disp(omega(z))
    
    %Time parmeters
    T=2*pi/omega(z);
    t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
    dt=t(2)-t(1); % Time resolution of the drive
    
    W_0=-1;  % Initial 0 winding number at omega(1)
    
    % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    clear singularities_0 q0 idx kx_aux ky_aux kx_idx ky_idx
       
    for j=1:length(lambda)
       
        % Time evolution
        for i=1:length(kx_gamma_K);
            
            U=Time_evolution_bulk_lambda_2T(delta,kx_gamma_K(i),ky_gamma_K(i),n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_Gamma_K_1(i)=HF_diag(1)/J; % Quasi energy bands
            eps_Gamma_K_2(i)=HF_diag(2)/J; % Quasi energy bands
            
        end
        for i=1:length(kx_KK);
            
            U=Time_evolution_bulk_lambda_2T(delta, kx_KK(i), ky_KK(i), n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_K_1(i)= HF_diag(1)/J; % Quasi energy bands
            eps_K_K_2(i)= HF_diag(2)/J; % Quasi energy bands
            
            
        end
        for i=1:length(kx_K_gamma);
            
            U=Time_evolution_bulk_lambda_2T(delta, kx_K_gamma(i), ky_K_gamma(i), n_t, dt,J_1,J_2,J_3, lambda(j)); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_Gamma_1(i)= HF_diag(1)/J; % Quasi energy bands
            eps_K_Gamma_2(i)= HF_diag(2)/J; % Quasi energy bands
            
        end
        
        % Gap calculation
        gap_0=min([abs(eps_Gamma_K_1),abs(eps_K_K_1),abs(eps_K_Gamma_1)]);
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
        
    end
    
    % Singularities and charges    
    if cont_0>1
        q0=top_charge_0_lambda(cont_0, singularities_0, step, n_t, delta,J_1,J_2,J_3); % Assign a charge to each singularity
        
    else
        q0=0;
    end   

    % Phase diagram
    for j=1:length(lambda)
        
        for i=1:length(q0(:,1))
            if lambda(j)==q0(i,1)
                W_0=W_0-q0(i,2); % Adding the topological charge at each singularity
            end
        end
        
        phase_diag(z,j,:)=[omega(z), lambda(j), W_0];
        
    end
    
end

%% Figures

figure(2)
close all
hold on
box on

for z=2:2:length(omega)-1
    
    for j=2:5:length(lambda)-1
        
        if phase_diag(z,j,3)==0 ;
            plot( lambda(j), omega(z), '.b','MarkerSize',8)
            
        else if phase_diag(z,j,3)==-1 ;
                plot(lambda(j), omega(z), '.c','MarkerSize',8)
                
            else if phase_diag(z,j,3)==1;
                    plot(lambda(j), omega(z), '.g','MarkerSize',8)
                    
                    
                else if phase_diag(z,j,3)==-2;
                        plot(lambda(j), omega(z), '.m','MarkerSize',8)
                        
                    else
                        plot(lambda(j), omega(z), '.k','MarkerSize',8)
                    end
                end
            end
        end
    end
end



%title(['Phase Diagram for \Delta/J=', num2str(delta)])
ax = gca;
ax.FontSize =16; 
xlabel('$\lambda$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\omega/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(2.5,6,'$[0,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(4.2,5.2,'$[0,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(2.9,5,'$[-1,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(4,6.7,'$[1,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')










