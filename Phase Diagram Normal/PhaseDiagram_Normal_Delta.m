%% PHASE DIAGRAM DELTA DRIVING

clear all
home
addpath('Functions')

% IMPORTANT COMMENTS: The drive is in delta, from high to low delta (low to high frequency) so
% charges must be subtracted for the winding numbers. Set the tolerance at
% around 0.1 when looking at anisotropic cases, when looking at isotropic
% cases 0.01 should be sufficient

%% Parameters
% Drive Family Physical Parameters
step=0.01; % Resolution of the frequency vector
omega=[1.5]; %  Frequency vector (units of J)
periods=2*pi./omega; % Period vector
delta=[1.1:-step:0]; % Potential offset vector (Sweeping delta in inverse sense, charges must be substracted not added!!)
J=1;
J_1=0.5;
J_2=1.5;
J_3=1;

% Computational parameters
tol=0.5; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vector used in the time evolution

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

%High Symmetry Points
% This part should be commented out if looking for the anisotropic case!!!
kx_gamma_K=0;   % From the Gamma point to K
ky_gamma_K=0;
kx_KK=l;        % From K to K'
ky_KK=0;
kx_K_gamma=s;   % From K to the Gamma point
ky_K_gamma=f;
momentumx=[kx_gamma_K, kx_KK, kx_K_gamma];     
momentumy=[ky_gamma_K, ky_KK, ky_K_gamma]; 
if J_1~=J || J_2~=J || J_3~=J
    disp('Anisotropic case only considering high symmetry points...')
end

% % High Symmetry Lines
% kx_K2_gamma=[s:-0.05:0];   % From K2 to the Gamma point
% ky_K2_gamma=tan(60*pi/180)*kx_K2_gamma;
% kx_gamma_K1=[0:0.05:l];    % From the Gamma point to K1
% ky_gamma_K1=0*kx_gamma_K1;
% kx_K1_K2=[l:-0.05:s];      % From K1 to K2
% ky_K1_K2=tan(120*pi/180)*kx_K1_K2-(l*tan(120*pi/180));
% kx_K2_K3=[s:-0.05:-s];     % From K2 to K3
% ky_K2_K3=f*ones(1,length(kx_K2_K3));
% kx_K3_K4=[-s:-0.05:-l];    % From K3 to K4
% ky_K3_K4=tan(60*pi/180)*kx_K3_K4+(l*tan(60*pi/180));
% kx_K4_gamma=[-l:0.05:0];   % From K4 to the Gamma point
% ky_K4_gamma=0*kx_K4_gamma;
% kx_gamma_K3=[0:-0.05:-s];  % From the Gamma point to K3
% ky_gamma_K3=-tan(60*pi/180)*kx_gamma_K3;
% kx_K5_gamma=[-s:0.05:0];   % From K5 to the Gamma point
% ky_K5_gamma=tan(60*pi/180)*kx_K5_gamma;
% kx_gamma_K6=[0:0.05:s];  % From the Gamma point to K6
% ky_gamma_K6=-tan(60*pi/180)*kx_gamma_K6;
% momentumx=[kx_K2_gamma, kx_gamma_K1, kx_K1_K2, kx_K2_K3, kx_K3_K4, kx_K4_gamma, kx_gamma_K3, kx_K5_gamma, kx_gamma_K6];     
% momentumy=[ky_K2_gamma, ky_gamma_K1, ky_K1_K2, ky_K2_K3, ky_K3_K4, ky_K4_gamma, ky_gamma_K3, ky_K5_gamma, ky_gamma_K6];  

%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.

for z=1:length(omega)
    
    disp(omega(z)) % Check

    % Time parameters
    T=periods(z); % Period for the j-th drive
    t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
    dt=t(2)-t(1); % Time resolution of the drive
       
    % Initial winding numbers
    W_0=0;  % Initial 0 winding number at omega(1)
    W_pi=0; % Initial pi winding number at omega(2)
    
   % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    cont_pi=1; % Number of times a singularity pi is detected  
    clear singularities_0 singularities_pi q0 qpi idx kx_idx ky_idx bands_1 bands_2
    
    for j=1:length(delta)
        
        % Time evolution
        for i=1:length(momentumx);
            
            U=Time_evolution_bulk(delta(j), momentumx(i), momentumy(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            bands_1(i)=real(HF_diag(1)); % Quasi energy bands
            bands_2(i)=real(HF_diag(2)); % Quasi energy bands
            
        end
        
        % Gap calculation
        gap_0=min(abs(bands_1));
        gap_pi=min(pi*ones(1,length(bands_1))-abs(bands_1));
        if  gap_0<tol
            idx = find(abs(bands_1)==gap_0); % Position of the gap
            kx_idx=momentumx(idx(1)); % x-momentum of the gap, idx(1) in case there are 2 simultaneous gaps
            ky_idx=momentumy(idx(1)); % y-momentum of the gap
            
            % Point in 3D parameter space where the band touching occurs
            singularities_0(cont_0,:) = [delta(j), kx_idx, ky_idx, dt, gap_0]; % We save dt for the S matrix
            cont_0=cont_0+1; % 1 extra singularity has been detected
        end
        if  gap_pi<tol
            
            idx = find(pi*ones(1,length(bands_1))-abs(bands_1)==gap_pi);
            kx_idx=momentumx(idx(1)); % x-momentum of the gap
            ky_idx=momentumy(idx(1)); % y-momentum of the gap
            h_0=-(abs(bands_1(idx(1)))+abs(bands_2(idx(1))))/2; % Calculation of h_0 in order to decompose on Pauli matrices
            
            
            % Point in 3D parameter space where the band touching occurs
            singularities_pi(cont_pi,:) = [delta(j), kx_idx, ky_idx, dt, h_0]; % We save dt for the S matrix
            cont_pi=cont_pi+1; % 1 extra singularity has been detected
        end
    end
     
    % Singularities and charges    
    if cont_0>1
        [q0, kx0, ky0]=top_charge_0_delta(cont_0, singularities_0, step, n_t, omega(z),J_1,J_2,J_3); % Assign a charge to each singularity
    else
        q0=[0,0];
    end  
    if cont_pi>1
        [qpi, kxpi, kypi]=top_charge_pi_delta(cont_pi, singularities_pi, step, n_t, omega(z),J_1,J_2,J_3); % Assign a charge to each singularity
    else
        qpi=[0,0];
    end  

    % Phase diagram
    for j=1:length(delta)
        
        for i=1:length(q0(:,1))
            if delta(j)==q0(i,1)
                W_0=W_0-q0(i,2); % Adding the topological charge at each singularity
            end
        end
        for i=1:length(qpi(:,1))
            if delta(j)==qpi(i,1)
                W_pi=W_pi-qpi(i,2);  % Adding the topological charge at each singularity
            end
        end
        
        phase_diag(j,z,:)=[delta(j), omega(z), W_0, W_pi];
        
    end
    
end

%% Figures

%load('phase_diag_data2.mat','phase_diag')   
%save('phase_diag_data2.mat','phase_diag')

figure(2)
close all
hold on
box on

for z=1:length(omega)
    
    for j=1:length(delta)
        
        if phase_diag(j,z,3)==0 && phase_diag(j,z,4)==0 ;
            plot( omega(z), delta(j), '.b','MarkerSize',5)
            
        else if phase_diag(j,z,3)==1 && phase_diag(j,z,4)==0 ;
                plot(omega(z), delta(j), '.g','MarkerSize',5)
                
            else if phase_diag(j,z,3)==0 && phase_diag(j,z,4)==1 ;
                    plot(omega(z), delta(j), '.c','MarkerSize',5)
                    
                    
                else if phase_diag(j,z,3)==1 && phase_diag(j,z,4)==1 ;
                        plot(omega(z), delta(j), '.r','MarkerSize',5)
                        
                    else if phase_diag(j,z,3)==1 && phase_diag(j,z,4)==-1 || phase_diag(j,z,3)==-1 && phase_diag(j,z,4)==1;
                            plot(omega(z), delta(j), '.k','MarkerSize',5)
                            
                        else if phase_diag(j,z,3)==-1 && phase_diag(j,z,4)==0 || phase_diag(j,z,3)==0 && phase_diag(j,z,4)==-1;
                                plot(omega(z), delta(j), '.y','MarkerSize',5)
                                                               
                            else
                                plot(omega(z), delta(j), '.m','MarkerSize',5)
                            end
                        end
                    end
                end
            end
        end
    end
end
title('Phase Diagram')
xlabel('$\omega/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\Delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')

ylim([0 2])
% xlim([0.7 3])



