%% PHASE DIAGRAM OMEGA DRIVE

clear all
home
addpath('Functions')

% IMPORTANT COMMENTS: The drive is in omega, from high to low frequency so
% charges must be summed up for the winding numbers. Set the tolerance at
% around 0.1 when looking at anisotropic cases, when looking at isotropic
% cases 0.01 should be sufficient


%% Parameters

% Drive Family Physical Parameters
step=0.001; % Resolution of the frequency vector
omega=[3:-step:1]; %  Frequency vector (units of J)
periods=2*pi./omega; % Period vector
delta=[0:0.1:2]; % Potential offset vector
J=1; % Every quantity in terms of J
J_1=1*J;
J_2=1*J;
J_3=1*J;

% Analysis of where the singularities occur in K-Space
cont_mom_0=1;
cont_mom_pi=1;
vec_mom_0=zeros(length(delta),4,3); % Storing charge and point in k-space
vec_mom_pi=zeros(length(delta),4,3); % Storing charge and point in k-space

% Computational parameters
tol=0.01; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vehomctor used in the time evolution

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


%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.

for z=1:length(delta)
    
    disp(delta(z)) % Check
    
    % Winding numbers at high frequency
    if abs(delta(z))>0.42
        W_0=0;  % Initial 0 winding number at omega(1)
    else
        W_0=1;  % Initial 0 winding number at omega(1)
    end
    
    
    % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    clear singularities_0 q0 idx kx_idx ky_idx bands_1 bands_2
    
    for j=1:length(omega)
        
        % Time parameters
        T=periods(j); % Period for the j-th drive
        t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,3*T/3,n_t)]; % Time vector of the drive
        dt=t(2)-t(1); % Time resolution of the drive
        
        % Time evolution
        for i=1:length(momentumx);
            
            U=Time_evolution_bulk_2T(delta(z), momentumx(i), momentumy(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
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
            singularities_0(cont_0,:) = [omega(j), kx_idx, ky_idx, dt, gap_0]; % We save dt for the S matrix
            cont_0=cont_0+1; % 1 extra singularity has been detected
        end
        
    end
    
    % Singularities and charges
    if cont_0>1
        [q0, kx0, ky0]=top_charge_0(cont_0, singularities_0, step, n_t, delta(z),J_1,J_2,J_3); % Assign a charge to each singularity
        vec_mom_0(cont_mom_0,1:length(q0(:,2)),:)=[q0(:,2), kx0, ky0]; % Map of the singulrities in the BZ
        cont_mom_0=cont_mom_0+1;
    else
        q0=[0];
    end
    
    % Phase diagram
    for j=1:length(omega)
        
        for i=1:length(q0(:,1))
            if omega(j)==q0(i,1)
                W_0=W_0+q0(i,2); % Adding the topological charge at each singularity
            end
        end
        
        phase_diag(z,j,:)=[delta(z), omega(j), W_0];
    end
    
   
    
end
    

%% Figures


figure(2)
close all
hold on
box on

for z=1:length(delta)
    
    for j=2:5:length(omega)-1
        
        if phase_diag(z,j,3)==0;
            plot( omega(j), delta(z), '.b','MarkerSize',8)
            
        else if phase_diag(z,j,3)==-1  ;
                plot(omega(j), delta(z), '.c','MarkerSize',8)              
                    
                else if phase_diag(z,j,3)==1  ;
                        plot(omega(j), delta(z), '.g','MarkerSize',8)
                                                  
                        else if phase_diag(z,j,3)==2;
                                plot(omega(j), delta(z), '.r','MarkerSize',8)
                                                               
                            else
                                plot(omega(j), delta(z), '.m','MarkerSize',8)
                            end
                        end
                    end
                end
            end
end
%title(['Phase Diagram for J_1=', num2str(J_1),' J_2=', num2str(J_2), ' J_3=',num2str(J_3)])
%legend('W_0','W_{\pi}','Location','Best')
xlabel('$\omega/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\Delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')

text(2.4,0.25,'$[1,0]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(2.7,1,'$[0,0]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(2,1.3,'$[0,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(1.5,0.5,'$[1,1]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(1.3,1.4,'$[1,0]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(1.05,1.8,'$[0,0]$','FontSize',25,'FontWeight','bold','Interpret','latex')
text(1.7,1.8,'$[0,0]$','FontSize',25,'FontWeight','bold','Interpret','latex')
ax = gca;
ax.FontSize =16; 
ax.XTick = [1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3];
ax.XTickLabel = {'1','1.2','1.4','1.6','1.8','2','2.2','2.4','2.6','2.8','3'};
ax.YTick = [0,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2];
ax.YTickLabel = {'0','0.2','0.4','0.6','0.8','1','1.2','1.4','1.6','1.8','2'};
%print -djpeg -r600 phasedigram1.png
% print -depsc phasediagram1.eps
