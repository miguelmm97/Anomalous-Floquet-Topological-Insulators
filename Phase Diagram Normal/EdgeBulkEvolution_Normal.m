% ANIMATED PLOT OF EDGE VS BULK STATES

% Here we make an animation of how the bulk states in the high symmetry
% lines and the edge states in the armchair system vary over the range of
% frequency, in order to be able to see the phase transitions more clearly

clear all
close all
addpath('Functions')
home

%% Parameters
% Drive Family Physical Parameters
step=0.1; % Resolution of the frequency vector
omega=[3:-step:1]; %  Frequency vector (units of J)
periods=2*pi./omega; % Period vector
delta_vec=[0]; % Potential offset 
J_1=1; % Hoping amplitudes (units of J)
J_2=1;   % Hoping amplitudes (units of J)
J_3=1;   % Hoping amplitudes (units of J)

% Computational parameters
n_t=100;

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% Edge Hamiltonian parameters
N_layers=32;  % Number of layers of the edge system (Armchair termination)
main_diag=ones(N_layers,1); % Parameters to construct the edge hamiltonian
sup_diag=ones(N_layers-1,1); % Parameters to construct the edge hamiltonian
inf_diag=ones(N_layers-1,1); % Parameters to construct the edge hamiltonian

% High Symmetry Lines
kx_K2_gamma=[s:-0.05:0];   % From K2 to the Gamma point
ky_K2_gamma=tan(60*pi/180)*kx_K2_gamma;
kx_gamma_K1=[0:0.05:l];    % From the Gamma point to K1
ky_gamma_K1=0*kx_gamma_K1;
kx_K1_K2=[l:-0.05:s];      % From K1 to K2
ky_K1_K2=tan(120*pi/180)*kx_K1_K2-(l*tan(120*pi/180));
kx_K2_K3=[s:-0.05:-s];     % From K2 to K3
ky_K2_K3=f*ones(1,length(kx_K2_K3));
kx_K3_K4=[-s:-0.05:-l];    % From K3 to K4
ky_K3_K4=tan(60*pi/180)*kx_K3_K4+(l*tan(60*pi/180));
kx_K4_gamma=[-l:0.05:0];   % From K4 to the Gamma point
ky_K4_gamma=0*kx_K4_gamma;
kx_gamma_K3=[0:-0.05:-s];  % From the Gamma point to K3
ky_gamma_K3=-tan(60*pi/180)*kx_gamma_K3;
kx_K5_gamma=[-s:0.05:0];   % From K5 to the Gamma point
ky_K5_gamma=tan(60*pi/180)*kx_K5_gamma;
kx_gamma_K6=[0:0.05:s];  % From the Gamma point to K6
ky_gamma_K6=-tan(60*pi/180)*kx_gamma_K6;

% Momentum range for edge states
kx_min=linspace(-pi,0,50); % Parallel momentum in the edge system
kx_max=linspace(0,pi,50); % Parallel momentum in the edge system
kx=[kx_min,kx_max];  % Parallel momentum in the edge system

% Momentum range for the full BZ
kx_BZ=[-l:0.1:l]; % Span the full BZ (really more than it as it is a square lxl)
ky_BZ=[-l:0.1:l];
K_BZ=[0,l];
K2_BZ=[f,s];

% Plotting parameters
p1=sqrt(s^2+f^2); % Distance covered in the brillouin zone to the end of  each line
p2=l;
p3=sqrt(f^2+(l-s)^2);
p4=2*s;
p5=p3;
p6=p2;
p7=p3;
p8=p3;
p9=p3;
k1=linspace(0, p1, length(kx_K2_gamma)); % Auxiliary vectors to plot the bad structure
k2=linspace(p1, p1+p2, length(kx_gamma_K1));
k3=linspace(p1+p2, p1+p2+p3, length(kx_K1_K2)); 
k4=linspace(p1+p2+p3, p1+p2+p3+p4, length(kx_K2_K3)); 
k5=linspace(p1+p2+p3+p4, p1+p2+p3+p4+p5, length(kx_K3_K4)); 
k6=linspace(p1+p2+p3+p4+p5, p1+p2+p3+p4+p5+p6, length(kx_K4_gamma)); 
k7=linspace(p1+p2+p3+p4+p5+p6, p1+p2+p3+p4+p5+p6+p7, length(kx_gamma_K3)); 
k8=linspace(p1+p2+p3+p4+p5+p6+p7, p1+p2+p3+p4+p5+p6+p7+p8, length(kx_K5_gamma)); 
k9=linspace(p1+p2+p3+p4+p5+p6+p7+p8, p1+p2+p3+p4+p5+p6+p7+p8+p9, length(kx_gamma_K6)); 

%% Animated plot comparing bulk and edge states
close all
for z=1:length(delta_vec)    
    for j=1:length(omega)
        
        % Potential offset parameters
        delta=delta_vec(z);
        diagonal=(delta/2)*ones(N_layers,1); % Construction edge hamiltonian
        
        % Time parameters
        T=periods(j);
        t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,3*T/3,n_t)]; % Time vector of the drive
        %t=[linspace(0,T/4,n_t),linspace(T/4,2*T/4,n_t),linspace(2*T/4,3*T/4,n_t),linspace(4*T/4,T,n_t)]; % Time vector of the drive
        dt=t(2)-t(1);
        
        % Edge states
        for i=1:length(kx);
            
            % Hamiltonian construction (Check hamiltonian edge in the notes)
            for s=1:N_layers
                if mod(s,2)==0
                    main_diag(s)=exp(1i*kx(i));
                end
            end
            gamma1=J_1*diag(main_diag,0);
            gamma2=J_2*diag(inf_diag,-1);
            gamma3=J_3*diag(sup_diag,1);
            
            H1=[diag(diagonal),gamma1;
                conj(transpose(gamma1)), diag(-diagonal)];
            H2=[diag(diagonal),gamma2;
                conj(transpose(gamma2)),diag(-diagonal)];
            H3=[diag(diagonal),gamma3;
                conj(transpose(gamma3)),diag(-diagonal)];
            
            % Time evolution operator and quasi-energy
            U1=(expm(-1i*dt*H1))^n_t;
            U2=(expm(-1i*dt*H2))^n_t;
            U3=(expm(-1i*dt*H3))^n_t;
            prod=U1*U2*U3;
            U_diag=eig(prod);
            
            % Quasi energy band structure
            for f=1:2*N_layers
                eps_edge(f,i)=1i*log(U_diag(f))./T;
            end
        end        
        % Bulk States
        for i=1:length(kx_K2_gamma);
            
            U=Time_evolution_bulk(delta, kx_K2_gamma(i),ky_K2_gamma(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K2_gamma_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K2_gamma_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_gamma_K1);
            
            U=Time_evolution_bulk(delta, kx_gamma_K1(i),ky_gamma_K1(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_gamma_K1_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_gamma_K1_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_K1_K2);
            
            U=Time_evolution_bulk(delta, kx_K1_K2(i),ky_K1_K2(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K1_K2_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K1_K2_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_K2_K3);
            
            U=Time_evolution_bulk(delta, kx_K2_K3(i),ky_K2_K3(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K2_K3_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K2_K3_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_K3_K4);
            
            U=Time_evolution_bulk(delta, kx_K3_K4(i),ky_K3_K4(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K3_K4_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K3_K4_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_K4_gamma);
            
            U=Time_evolution_bulk(delta, kx_K4_gamma(i),ky_K4_gamma(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K4_gamma_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K4_gamma_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_gamma_K3);
            
            U=Time_evolution_bulk(delta, kx_gamma_K3(i),ky_gamma_K3(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_gamma_K3_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_gamma_K3_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_K5_gamma);
            
            U=Time_evolution_bulk(delta, kx_K5_gamma(i),ky_K5_gamma(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_K5_gamma_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_K5_gamma_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end
        for i=1:length(kx_gamma_K6);
            
            U=Time_evolution_bulk(delta, kx_gamma_K6(i),ky_gamma_K6(i), n_t, dt, J_1, J_2, J_3); % Floquet Hamiltonian (H*T)
            HF=1i*logm(U);
            HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
            eps_gamma_K6_1(i)=real(HF_diag(1))./T; % Quasi energy bands
            eps_gamma_K6_2(i)=real(HF_diag(2))./T; % Quasi energy bands
            
        end


        
        figure (1)
        suptitle(['  \omega=',num2str(2*pi/T), ' \Delta /J =', num2str(delta)])
        subplot(1,2,1)
        box on
        hold on
        for j=1:2*N_layers
            plot(kx,eps_edge(j,:)*T,'.b','MarkerSize',3)
        end   % Plot all the edge states
        hold off
        title('Edge Spectrum','FontSize',12,'FontWeight','bold')
        xlabel('$k_{\parallel}a$','FontSize',18,'FontWeight','bold','Interpret','latex')
        ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
        ax=gca;
        ax = gca;
        ax.FontSize = 12;
        xlim([-pi pi])
        ylim([-pi pi])
        
        subplot(1,2,2)
        box on
        hold on
        plot(k1,eps_K2_gamma_1*T,'b.',k1,eps_K2_gamma_2*T,'b.')
        plot(k2,eps_gamma_K1_1*T,'.b',k2,eps_gamma_K1_2*T,'.b')
        plot(k3, eps_K1_K2_1*T,'.b',k3, eps_K1_K2_2*T, '.b')
        plot(k4, eps_K2_K3_1*T,'.b',k4, eps_K2_K3_2*T, '.b')
        plot(k5, eps_K3_K4_1*T,'.b',k5, eps_K3_K4_2*T, '.b')
        plot(k6, eps_K4_gamma_1*T,'.b',k6, eps_K4_gamma_2*T, '.b')
        plot(k7, eps_gamma_K3_1*T,'.b',k7, eps_gamma_K3_2*T, '.b')
        plot(k8, eps_K5_gamma_1*T,'.b',k8, eps_K5_gamma_2*T, '.b')
        plot(k9, eps_gamma_K6_1*T,'.b',k9, eps_gamma_K6_2*T, '.b')
        
        title('Bulk Spectrum','FontSize',12,'FontWeight','bold')
        xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
        ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
        hold off
        ax=gca;
        ax.XTick=[0,p1, p1+p2, p1+p2+p3, p1+p2+p3+p4, p1+p2+p3+p4+p5, p1+p2+p3+p4+p5+p6, p1+p2+p3+p4+p5+p6+p7, p1+p2+p3+p4+p5+p6+p7+p8, p1+p2+p3+p4+p5+p6+p7+p8+p9 ];
        ax.XTickLabel={'K_2','\Gamma','K_1','K_2','K_3','K_4','\Gamma','K3','\Gamma','K_6'}
        ax = gca;
        ax.FontSize = 12;
        xlim([0 p1+p2+p3+p4+p5+p6+p7+p8+p9 ])
        ylim([-pi pi])
        
        set(gcf,'Position',[900 600 1000 400]) % Sets position and width
        drawnow  % Animation
        clf      % Clear after each step in the loop
             
    end   
end

%% Animated plot of bulk evolution in the full BZ
close all
for z=1:length(delta_vec)
    
    delta=delta_vec(z);
    
    for w=1:length(omega)
        
        T=periods(w);
        t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
        dt=t(2)-t(1);
        
        for i=1:length(kx_BZ);
            for j=1:length(ky_BZ)
                                
                U=Time_evolution_bulk(delta,kx_BZ(i),ky_BZ(j),n_t, dt, J_1,J_2,J_3); % Floquet Hamiltonian (H*T)
                HF=1i*logm(U);
                HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                eps_BZ_1(i,j)=real(HF_diag(1)); % Quasi energy bands
                eps_BZ_2(i,j)=real(HF_diag(2)); % Quasi energy bands
                
            end
        end
        
        
        figure(1)               
        suptitle(['  \omega=',num2str(2*pi/T), ' \Delta /J =', num2str(delta)])
        subplot(1,2,2)
        hold on
        box on
        contour(kx_BZ,ky_BZ, abs(eps_BZ_1),'ShowText','on');
        scatter(0,0,'r','filled')
        scatter(K_BZ(1),K_BZ(2),'r','filled')
        scatter(K2_BZ(1),K2_BZ(2),'r','filled')
        scatter(-K_BZ(1),K_BZ(2),'r','filled')
        scatter(-K2_BZ(1),K2_BZ(2),'r','filled')
        scatter(K_BZ(1),-K_BZ(2),'r','filled')
        scatter(K2_BZ(1),-K2_BZ(2),'r','filled')
        scatter(-K2_BZ(1),-K2_BZ(2),'r','filled')
        caxis([0 pi])
        xlim([-f f])
        ylim([-l l])
        xlabel('$k_x a$','FontSize',18,'FontWeight','bold','Interpret','latex')
        ylabel('$k_y a$','FontSize',18,'FontWeight','bold','Interpret','latex')
        hold off
        
        subplot(1,2,1)
        hold on
        surf(kx_BZ,ky_BZ,abs(eps_BZ_1));
        surf(kx_BZ,ky_BZ,-abs(eps_BZ_2));
        box on
        view([15 10]) % Sets the view of the camera sight
        caxis([-pi pi])
        xlim([-f f])
        ylim([-l l])
        zlim([-pi pi])
        xlabel('$k_x a$','FontSize',18,'FontWeight','bold','Interpret','latex')
        ylabel('$k_y a$','FontSize',18,'FontWeight','bold','Interpret','latex')
        zlabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
        hold off
        
        set(gcf,'Position',[500 600 1300 400])
        
        %set(gcf,'Position',[900 600 1000 400]) % Sets position and width
        drawnow  % Animation
        clf      % Clear after each step in the loop
        
    end
end


