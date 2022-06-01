%% PHASE DIAGRAM FULL BZ OMEGA DRIVE

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
omega=[4.5:-step:1]; %  Frequency vector (units of J)
periods=2*pi./omega; % Period vector
delta=[0:0.3:2]; % Potential offset vector
J=1; % Every quantity in terms of J
J_1=1*J;
J_2=1*J;
J_3=3*J;

% Analysis of where the singularities occur in K-Space
cont_mom_0=1;
cont_mom_pi=1;
vec_mom_0=zeros(length(delta),4,3); % Storing charge and point in k-space
vec_mom_pi=zeros(length(delta),4,3); % Storing charge and point in k-space

% Computational parameters
tol=0.1; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vehomctor used in the time evolution

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% Full BZ
momentumx=[-l:0.1:l]; % Spans a bit more than the full BZ but it can be easily distinguished.
momentumy=[-l:0.1:l];


%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.

for z=1:length(delta)
    
    disp(delta(z)) % Check
    
    % Winding numbers at high frequency    
    if delta(z)>0.45
        W_0=0;  % Initial 0 winding number at omega(1)
        W_pi=0; % Initial pi winding number at omega(2)
    else
        W_0=1;
        W_pi=0;
    end
       
    
    % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    cont_pi=1; % Number of times a singularity pi is detected  
    clear singularities_0 singularities_pi q0 qpi idx1 idx2 kx_idx ky_idx bands_1 bands_2
          
    for j=1:length(omega)    
      
        % Time parameters
        T=periods(j); % Period for the j-th drive
        t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,3*T/3,n_t)]; % Time vector of the drive
        dt=t(2)-t(1); % Time resolution of the drive
        
        % Time evolution
        for i=1:length(momentumx);
            for s=1:length(momentumy)
            U=Time_evolution_bulk(delta(z), momentumx(i), momentumy(s), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            bands_1(i,s)=real(HF_diag(1)); % Quasi energy bands
            bands_2(i,s)=real(HF_diag(2)); % Quasi energy bands
            
            end
        end
             
        % Gap calculation                                                                                   
        gap_0=min(min(abs(bands_1)));
        gap_pi=min(min((pi*ones(length(bands_1(:,1)))-abs(bands_1))));               
        if  gap_0<tol
            [idx1, idx2] = find(abs(bands_1)==gap_0); % Position of the gap
            kx_idx=momentumx(idx1(1)); % x-momentum of the gap, idx(1) in case there are 2 simultaneous gaps
            ky_idx=momentumy(idx2(1)); % y-momentum of the gap
            
            % Point in 3D parameter space where the band touching occurs
            singularities_0(cont_0,:) = [omega(j), kx_idx, ky_idx, dt, gap_0]; % We save dt for the S matrix
            cont_0=cont_0+1; % 1 extra singularity has been detected
        end   
        if  gap_pi<tol
            
            [idx1, idx2] = find((pi*ones(length(bands_1(:,1)))-abs(bands_1))==gap_pi);
            kx_idx=momentumx(idx1(1)); % x-momentum of the gap
            ky_idx=momentumy(idx2(1)); % y-momentum of the gap            
            h_0=-(abs(bands_1(idx1(1),idx2(1)))+abs(bands_2(idx1(1),idx2(1))))/2; % Calculation of h_0 in order to decompose on Pauli matrices
            
            
            % Point in 3D parameter space where the band touching occurs
            singularities_pi(cont_pi,:) = [omega(j), kx_idx, ky_idx, dt, h_0]; % We save dt for the S matrix
            cont_pi=cont_pi+1; % 1 extra singularity has been detected
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
    if cont_pi>1
        [qpi, kxpi, kypi]=top_charge_pi(cont_pi, singularities_pi, step, n_t, delta(z),J_1,J_2,J_3); % Assign a charge to each singularity
        vec_mom_pi(cont_mom_pi,1:length(qpi(:,2)),:)=[qpi(:,2), kxpi, kypi];  % Map of the singulrities in the BZ
        cont_mom_pi=cont_mom_pi+1;
    else
        qpi=[0];
    end  

    % Phase diagram
    for j=1:length(omega)
        
        for i=1:length(q0(:,1))
            if omega(j)==q0(i,1)
                W_0=W_0+q0(i,2); % Adding the topological charge at each singularity
            end
        end
        for i=1:length(qpi(:,1))
            if omega(j)==qpi(i,1)
                W_pi=W_pi+qpi(i,2);  % Adding the topological charge at each singularity
            end
        end
        
        phase_diag(z,j,:)=[delta(z), omega(j), W_0, W_pi];
        
    end
    
end

%% Figures

figure(2)
close all
hold on
box on

for z=1:length(delta)
    
    for j=1:1:10:length(omega)
        
        if phase_diag(z,j,3)==0 && phase_diag(z,j,4)==0 ;
            plot( omega(j), delta(z), '.b','MarkerSize',5)
            
        else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==0 ;
                plot(omega(j), delta(z), '.g','MarkerSize',5)
                
            else if phase_diag(z,j,3)==0 && phase_diag(z,j,4)==1 ;
                    plot(omega(j), delta(z), '.c','MarkerSize',5)
                    
                    
                else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==1 ;
                        plot(omega(j), delta(z), '.r','MarkerSize',5)
                        
                    else if phase_diag(z,j,3)==1 && phase_diag(z,j,4)==-1 || phase_diag(z,j,3)==-1 && phase_diag(z,j,4)==1;
                            plot(omega(j), delta(z), '.k','MarkerSize',5)
                            
                        else if phase_diag(z,j,3)==-1 && phase_diag(z,j,4)==0 || phase_diag(z,j,3)==0 && phase_diag(z,j,4)==-1;
                                plot(omega(j), delta(z), '.y','MarkerSize',5)
                                                               
                            else
                                plot(omega(j), delta(z), '.m','MarkerSize',5)
                            end
                        end
                    end
                end
            end
        end
    end
end
title(['Phase Diagram for J_1=', num2str(J_1),' J_2=', num2str(J_2), ' J_3=',num2str(J_3)])
%legend('W_0','W_{\pi}','Location','Best')
xlabel('$\omega/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')

ylim([0 2])
% xlim([0.7 3])











