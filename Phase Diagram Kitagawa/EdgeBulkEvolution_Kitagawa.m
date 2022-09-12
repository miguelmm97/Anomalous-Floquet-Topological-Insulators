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
omega_vec=[7:-step:4.5]; %  Frequency (units of J)
periods=2*pi./omega_vec; % Period 
delta_vec=[2]; % Potential offset vector
lambda_vec=[1.3]; % Degree of anisotropy (lambda=1 identified with high frequency static regime)
J=1; J_1=J; J_2=J; J_3=J; % Every quantity in terms of J
% drive='lambda';
drive='omega';

% Computational parameters
tol=0.01; % Band gap tolerance (default 0.01)
n_t=1000; % Resolution of the time vector used in the time evolution

% Reciprocal Lattice Parameters (Check sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% Hamiltonian parameters edge
N_layers=30;
main_diag=ones(N_layers,1);
sup_diag=ones(N_layers-1,1);
inf_diag=ones(N_layers-1,1);

% Momentum range edge
kx_min=linspace(-pi,0,100);
kx_max=linspace(0,pi,100);
kx=[kx_min,kx_max];  

% High Symmetry Lines
kx_K2_gamma=[s:-0.01:0];   % From K2 to the Gamma point
ky_K2_gamma=tan(60*pi/180)*kx_K2_gamma;
kx_gamma_K1=[0:0.01:l];    % From the Gamma point to K1
ky_gamma_K1=0*kx_gamma_K1;
kx_K1_K2=[l:-0.01:s];      % From K1 to K2
ky_K1_K2=tan(120*pi/180)*kx_K1_K2-(l*tan(120*pi/180));

% Plotting parameters
p1=sqrt(s^2+f^2); % Distance covered in the brillouin zone to the end of  each line
p2=l;
p3=sqrt(f^2+(l-s)^2);
k1=linspace(0, p1, length(kx_K2_gamma)); % Auxiliary vectors to plot the bad structure
k2=linspace(p1, p1+p2, length(kx_gamma_K1));
k3=linspace(p1+p2, p1+p2+p3, length(kx_K1_K2)); 


%% Animated plot comparing bulk and edge states
close all

for z=1:length(delta_vec)
    
    switch drive
        case 'lambda'
            
            for j=1:length(lambda_vec)
                
                % Potential offset parameters
                delta=delta_vec(z);
                lambda=lambda_vec(j);
                T=periods;
                t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
                dt=t(2)-t(1); % Time resolution of the drive
                diagonal=(delta/2)*ones(N_layers,1); % Construction edge hamiltonian
                
                % Edge states
                for i=1:length(kx);
                    
                    % Hamiltonian construction
                    for s=1:N_layers
                        if mod(s,2)==0
                            main_diag(s)=exp(1i*kx(i));
                        end
                    end
                    gamma1=J_1*diag(main_diag,0);
                    gamma2=J_2*diag(inf_diag,-1);
                    gamma3=J_3*diag(sup_diag,1);
                    
                    H1=[diag(diagonal),(lambda*gamma1)+gamma2+gamma3;
                        conj(transpose((lambda*gamma1)+gamma2+gamma3)), diag(-diagonal)];
                    
                    H2=[diag(diagonal),gamma1+(lambda*gamma2)+gamma3;
                        conj(transpose(gamma1+(lambda*gamma2)+gamma3)),diag(-diagonal)];
                    
                    H3=[diag(diagonal),gamma1+gamma2+(lambda*gamma3);
                        conj(transpose(gamma1+gamma2+(lambda*gamma3))),diag(-diagonal)];
                    
                    
                    % Time evolution operator and quasi-energy
                    U1=(expm(-1i*dt*H1))^n_t;
                    U2=(expm(-1i*dt*H2))^n_t;
                    U3=(expm(-1i*dt*H3))^n_t;
                    prod=U1*U2*U3;
                    U_diag=eig(prod);
                    
                    for f=1:2*N_layers
                        eps_edge(f,i)=1i*log(U_diag(f));
                    end
                end
                % Bulk States
                for i=1:length(kx_K2_gamma);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_K2_gamma(i),ky_K2_gamma(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_K2_gamma_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_K2_gamma_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                for i=1:length(kx_gamma_K1);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_gamma_K1(i),ky_gamma_K1(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_gamma_K1_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_gamma_K1_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                for i=1:length(kx_K1_K2);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_K1_K2(i),ky_K1_K2(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_K1_K2_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_K1_K2_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                
                figure (1)
                suptitle(['  \lambda=',num2str(lambda), ' \Delta /J =', num2str(delta)])
                subplot(1,2,1)
                box on
                hold on
                for j=1:2*N_layers
                    plot(kx,eps_edge(j,:),'.b','MarkerSize',3)
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
                
                title('Bulk Spectrum','FontSize',12,'FontWeight','bold')
                xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
                ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
                hold off
                ax=gca;
                ax.XTick=[0,p1, p1+p2, p1+p2+p3];
                ax.XTickLabel={'K_2','\Gamma','K_1'}
                ax = gca;
                ax.FontSize = 12;
                xlim([0 p1+p2+p3 ])
                ylim([-pi pi])
                
                
                set(gcf,'Position',[900 600 1000 400]) % Sets position and width
                drawnow  % Animation
                clf      % Clear after each step in the loop
                
            end
            
        case 'omega'
            for j=1:length(omega_vec)
                              
                % Drive parameters
                lambda=lambda_vec;
                delta=delta_vec(z);
                T=periods(j);
                t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)]; % Time vector of the drive
                dt=t(2)-t(1); % Time resolution of the drive
                diagonal=(delta/2)*ones(N_layers,1); % Construction edge hamiltonian
                               
                % Edge states
                for i=1:length(kx);
                    
                    % Hamiltonian construction
                    for s=1:N_layers
                        if mod(s,2)==0
                            main_diag(s)=exp(1i*kx(i));
                        end
                    end
                    gamma1=J_1*diag(main_diag,0);
                    gamma2=J_2*diag(inf_diag,-1);
                    gamma3=J_3*diag(sup_diag,1);
                    
                    H1=[diag(diagonal),(lambda*gamma1)+gamma2+gamma3;
                        conj(transpose((lambda*gamma1)+gamma2+gamma3)), diag(-diagonal)];
                    
                    H2=[diag(diagonal),gamma1+(lambda*gamma2)+gamma3;
                        conj(transpose(gamma1+(lambda*gamma2)+gamma3)),diag(-diagonal)];
                    
                    H3=[diag(diagonal),gamma1+gamma2+(lambda*gamma3);
                        conj(transpose(gamma1+gamma2+(lambda*gamma3))),diag(-diagonal)];
                    
                    
                    % Time evolution operator and quasi-energy
                    U1=(expm(-1i*dt*H1))^n_t;
                    U2=(expm(-1i*dt*H2))^n_t;
                    U3=(expm(-1i*dt*H3))^n_t;
                    prod=U1*U2*U3;
                    U_diag=eig(prod);
                    
                    for f=1:2*N_layers
                        eps_edge(f,i)=1i*log(U_diag(f));
                    end
                end
                % Bulk States
                for i=1:length(kx_K2_gamma);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_K2_gamma(i),ky_K2_gamma(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_K2_gamma_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_K2_gamma_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                for i=1:length(kx_gamma_K1);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_gamma_K1(i),ky_gamma_K1(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_gamma_K1_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_gamma_K1_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                for i=1:length(kx_K1_K2);
                    
                    U=Time_evolution_bulk_lambda(delta, kx_K1_K2(i),ky_K1_K2(i), n_t, dt, J_1, J_2, J_3, lambda); % Floquet Hamiltonian (H*T)
                    HF=1i*logm(U);
                    HF_diag=eig(HF); % Diagonal Floquet Hamiltonian
                    eps_K1_K2_1(i)=real(HF_diag(1))./T; % Quasi energy bands
                    eps_K1_K2_2(i)=real(HF_diag(2))./T; % Quasi energy bands
                    
                end
                
                figure (1)
                suptitle(['  \omega/J=',num2str(omega_vec(j)), ' \Delta /J =', num2str(delta)])
                subplot(1,2,1)
                box on
                hold on
                for j=1:2*N_layers
                    plot(kx,eps_edge(j,:),'.b','MarkerSize',3)
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
                
                title('Bulk Spectrum','FontSize',12,'FontWeight','bold')
                xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
                ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
                hold off
                ax=gca;
                ax.XTick=[0,p1, p1+p2, p1+p2+p3];
                ax.XTickLabel={'K_2','\Gamma','K_1'}
                ax = gca;
                ax.FontSize = 12;
                xlim([0 p1+p2+p3 ])
                ylim([-pi pi])
                
                
                set(gcf,'Position',[900 600 1000 400]) % Sets position and width
                drawnow  % Animation
                clf      % Clear after each step in the loop
                
            end
            
    end
end


