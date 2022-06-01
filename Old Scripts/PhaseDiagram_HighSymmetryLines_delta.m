%% PHASE DIAGRAM FLOQUET GRAPHENE

clear all
home
addpath('Functions')

% Drive Family Physical Parameters
step=0.01; % Resolution of the frequency vector
omega=[0.52]; %  Frequency vector (units of J)
periods=2*pi./omega; % Period vector
delta=[1.6:-step:0]; % Potential offset vector (Sweeping delta in inverse sense, charges must be substracted not added!!)
J_1=1.1;
J_2=1;
J_3=1;

% Analysis of where the singularities occur in K-Space
% cont_mom_0=1;
% cont_mom_pi=1;
% vec_mom_0=zeros(length(delta),4,3); % Storing charge and point in k-space
% vec_mom_pi=zeros(length(delta),4,3); % Storing charge and point in k-space

% Computational parameters
tol=0.1; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vector used in the time evolution

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% High Symmetry Lines 
kx_gamma_K=[0:0.1:l];    % From the Gamma point to K
ky_gamma_K=0*kx_gamma_K;
kx_KK=[l:-0.1:s];        % From K to K'
ky_KK=tan(120*pi/180)*kx_KK-(l*tan(120*pi/180));
kx_K_gamma=[s:-0.1:0];   % From K to the Gamma point
ky_K_gamma=tan(60*pi/180)*kx_K_gamma;

%% Family of drives lambda=omega

% We select an offset potential and we sweep the quasi energy spectrum of
% the system from high-frequency trivial topology (0,0) to low frequency,
% searching for band touchings and calculating the associated topological
% charge with the "efficient" scheme.

for z=1:length(omega)
    
    omega(z) % Check

    % Time parameters
    T=periods(z); % Period for the j-th drive
    t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
    dt=t(2)-t(1); % Time resolution of the drive
    
    W_0=0;  % Initial 0 winding number at omega(1)
    W_pi=0; % Initial pi winding number at omega(2)  

    % Definition and clearing of variables involved in the loop
    cont_0=1; % Number of times a singularity 0 is detected
    cont_pi=1; % Number of times a singularity pi is detected
    singularities_0=[];
    singularities_pi=[];
    q0=[];
    qpi=[];
    idx=0;
    kx_aux=[]; 
    ky_aux=[];
    kx_idx=[]; 
    ky_idx=[];
       
    for j=1:length(delta)
        j
               
        % Time evolution
        for i=1:length(kx_gamma_K);
            
            U=Time_evolution_bulk(delta(j),kx_gamma_K(i),ky_gamma_K(i),n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_Gamma_K_1(i)=HF_diag(1); % Quasi energy bands
            eps_Gamma_K_2(i)=HF_diag(2); % Quasi energy bands
            
        end
        for i=1:length(kx_KK);
            
            U=Time_evolution_bulk(delta(j), kx_KK(i), ky_KK(i), n_t, dt,J_1,J_2,J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_K_1(i)= HF_diag(1); % Quasi energy bands
            eps_K_K_2(i)= HF_diag(2); % Quasi energy bands
            
            
        end
        for i=1:length(kx_K_gamma);
            
            U=Time_evolution_bulk(delta(j), kx_K_gamma(i), ky_K_gamma(i), n_t, dt,J_1,J_2,J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K_Gamma_1(i)= HF_diag(1); % Quasi energy bands
            eps_K_Gamma_2(i)= HF_diag(2); % Quasi energy bands
            
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
            singularities_0(cont_0,:) = [delta(j), kx_idx, ky_idx, dt, gap_0]; % We save dt for the S matrix
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
            singularities_pi(cont_pi,:) = [delta(j), kx_idx, ky_idx, dt, h_0]; % We save dt for the S matrix
            cont_pi=cont_pi+1; % 1 extra singularity has been detected
        end
        
    end
    
    % Singularities and charges
    
    if cont_0>1
        [q0, kx0, ky0]=top_charge_0_delta(cont_0, singularities_0, step, n_t, omega(z),J_1,J_2,J_3);
%         vec_mom_0(cont_mom_0,1:length(q0(:,2)),:)=[q0(:,2), kx0, ky0];
%         cont_mom_0=cont_mom_0+1;
    else
        q0=[0,0];
    end
    
    
    if cont_pi>1
        [qpi, kxpi, kypi]=top_charge_pi_delta(cont_pi, singularities_pi, step, n_t, omega(z),J_1,J_2,J_3);
%         vec_mom_pi(cont_mom_pi,1:length(qpi(:,2)),:)=[qpi(:,2), kxpi, kypi];
%         cont_mom_pi=cont_mom_pi+1;
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
%legend('W_0','W_{\pi}','Location','Best')
xlabel('$\omega/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$\delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(2.4,0.2,'$[1,0]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(2.5,1,'$[0,0]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(2,1.3,'$[0,1]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(1.5,0.5,'$[1,1]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(1.3,1.4,'$[1,0]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(1.1,1.6,'$[0,0]$','FontSize',18,'FontWeight','bold','Interpret','latex')
text(1.7,1.8,'$[0,0]$','FontSize',18,'FontWeight','bold','Interpret','latex')
%text(0.75,0.1,'$[0,1]$','FontSize',18,'FontWeight','bold','Interpret','latex')

ylim([0 2])
% xlim([0.7 3])


% figure(3)
% hold on
% box on
% 
% for z=1:length(delta)
%     plot3(0,0,delta(z),'.k','MarkerSize',10)
%     plot3(2.4184,0,delta(z),'.k','MarkerSize',10)
%     plot3(1.2092,2.094,delta(z),'.k','MarkerSize',10)
%     text(0.3,0,delta(z),'$\Gamma$','FontSize',18,'FontWeight','bold','Interpret','latex')
%     text(2.6,0,delta(z),'$K$','FontSize',18,'FontWeight','bold','Interpret','latex')
%     text(1.3,2.094,delta(z),'$K_2$','FontSize',18,'FontWeight','bold','Interpret','latex')
%     for j=1:2
%         if vec_mom_0(z,j,1)==1
%             plot3(vec_mom_0(z,j,2),vec_mom_0(z,j,3),delta(z),'^r','MarkerSize',10)
%         else if vec_mom_0(z,j,1)==-1
%             plot3(vec_mom_0(z,j,2),vec_mom_0(z,j,3),delta(z),'^b','MarkerSize',10)
%             end
%         end
%     end
% end
% for z=1:length(delta)
%     for j=1:2
%         if vec_mom_pi(z,j,1)==1
%             plot3(vec_mom_pi(z,j,2),vec_mom_pi(z,j,3),delta(z),'or','MarkerSize',15)
%         else if vec_mom_pi(z,j,1)==-1
%             plot3(vec_mom_pi(z,j,2),vec_mom_pi(z,j,3),delta(z),'ob','MarkerSize',15)
%             end
%         end
%     end
% end  
% title('Singularities in K-Space')
% %legend('W_0','W_{\pi}','Location','Best')
% xlabel('$k_x/a$','FontSize',18,'FontWeight','bold','Interpret','latex')
% ylabel('$k_y/a$','FontSize',18,'FontWeight','bold','Interpret','latex')
% zlabel('$\delta/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
% 
% xlim([-3 3])
% ylim([-3 3])

