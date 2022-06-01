% WAVE-PACKET DYNAMICS IN THE ANOMALOUS PHASE
% Here we develop the dynamics of wavepackets initialised in the edges of a
% finite geometry in any topological phase, in order to see how the different
% chiralities and nature of the two possible types of edge sttes influence
% the time evolution of the wave-packets. We can select Rudner or Kitagawa
% drives, strip or cyllinder geometries, to apply a flattening to the
% spectrum, to produce a movie of the evolution or not, and to introduce
% impurities in the simulated sample.
home
close all
clear all
addpath('functions')

%% Parameters

%%%%%%%%%%%%%%%%%%%%% Code parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

drive='Kitagawa';                                                             % 'Kitagawa' % Select driving approach
geometry= 'St';% 'Cyllinder' % 'Circle'                                  % Select geometry
evolution= 'Stroboscopic'; % 'Exact'                                        % Evolution for the Rudner drive
flattening= 'Off'; % Off                                                    % Apply flattening to the quasienergy spectrum
impurities= 'Off';                                                          % Put impurities in the sample
movie= 'Off';                                                               % If we want to produce a movie of the evolution

%%%%%%%%%%%%%%%%%%%% Phase parameter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega= 5;    % 6.5 7 5 (Sweet spot Rudner drive w=4/3)                    % Frequency of the drive ([J])
delta_plot=0.5;
delta=2; % 2 3                                                              % Potential offset ([J])
lambda=3.4;   % 2.5 4.5 3.4                                                 % Degree of anisotropy
J=1;                                                                        % Coupling constant
J_1=J; J_2=J; J_3=J;                                                        % Hoping amplitudes ([J])

%%%%%%%%%%%%%%%%%%%% Time parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_t=100;                                                                    % Time discretisation
T=2*pi/omega;                                                               % Period
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];      % Time grid
dt=t(2)-t(1);                                                               % Time step

%%%%%%%%%%%%%%%%% Number of layers of the geometry %%%%%%%%%%%%%%%%%%%%%%%%

Nx=32; % Necessarilly mod(Nx,4)=0 for zig-zag edge along y                  % Number of x sites
Ny=15; % Necessarilly odd for armchair edge along x                         % Number of y sites
Ns=Nx*Ny/2;   Ns_real=Ns;                                                   % Total number of states
y_sep=sqrt(3)/2;                                                            % Separation between y layers so that every site is equidistant from the others in units of a=1
x_sep=3;  layers_x=Nx*3/8;                                                  % x length of each cell
center=[Nx*3/8, y_sep*(Ny-1)/2];                                            % Center of the circle geometry
Radius=Nx*3/8;                                                              % Radius of the circle geometry

%%%%%%%%%%%%%%%%%%%%%%%% Wave-packet parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%

x0=((Nx/4)*x_sep/2); y0=y_sep*(Ny-1);                                              % Initial position of the wave-packet's center
sigma_x=3; sigma_y=1;                                                      % Widths
t_strob=60*T; t_exact=30*T;                                                 % Final time of the packet evolution

%%%%%%%%%%%%%%%%%%%%%%% Reciprocal Lattice Parameters %%%%%%%%%%%%%%%%%%%%%

l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

%%%%%%%%%%%%%%%%%%%%% High symmetry points in the BZ %%%%%%%%%%%%%%%%%%%%%%

gamma_BZ=[0,0];
K_BZ=[l,0];
K2_BZ=[f,s];

%%%%%%%%%%%%%%%%%%%%%%%%%% Declaration of variables %%%%%%%%%%%%%%%%%%%%%%%

position_matx=zeros(Ns,3);                                                   % 1st colum state number, 2nd x position, 3rd y position
Psi_0=zeros(Ns,1);                                                           % Initial wavepacket
overlap_0=zeros(Ns,1);                                                       % Initial overlap of the wave packet
total_prob_0=zeros(Ns,1);                                                    % Total prob for 0 gap
total_prob_pi=zeros(Ns,1);                                                   % Total prob pi gap
step=0;                                                                      % Step counter in the evolution of the Rudner drive
avoided_list=0;                                                              % Thrown away points in the circle geometry
HF_pi = zeros(Ns); HF_0=zeros(Ns);                                           % Definition for hamiltonians at both gaps for symmetry connection


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Declaration of functions %%%%%%%%%%%%%%%%%%%%%

gauss_weight=@(x, y, sigma_x, sigma_y)  ...                                  % Gaussian weight of the wave-packet
        (1/(2*pi*sigma_x*sigma_y))*exp(-((x-x0)^2)/(2*sigma_x^2))* exp(-((y-y0)^2)/(2*sigma_y^2));       
conv_factor=@(value, maximum) value/maximum;                                 % Conversion factor to plot probability densities
distance=@(vector) sqrt((vector(1)-center(1))^2+(vector(2)-center(2))^2);    % Norm function

%% State-position assignment and Hamiltonian Construction
% Here we make explicit the assignation of each state to its possition in
% the finite strip we are considering, and build the hamiltonian out of the
% hopping amplitudes, this is just as a legend to chech how we label each site 
% and each hoping amplitude. REMEMBER, POSITION IS ALWAYS IN UNITS OF a=1,
% where a is the lattice constant. From each point there is always 3a
% up to the same point in the next cell. If we want a cyllinder geometry we
% needn't change anything in this section, because this is just to plot the
% honeycomb grid and check the hamiltonian is fine, but after that, we
% won't plot the connection between the first and last layer.


% Construction of the position matrix (n_state, xpos, ypos)
position_matx=position_matrix(Nx, Ny);                                       % State position and labelling
limx=max(position_matx(:,2))+1;  limy=max(position_matx(:,3))+1;             % Usual limits of the strip to plot later on
H=Finite_Hamiltonian( delta_plot, 1, 2, 3, Nx, Ny, position_matx);           % Test hamiltonian strip
%H=Finite_Hamiltonian_Cyllinder( delta, 1, 2, 3, Nx, Ny, position_matx);     % Test hamiltonian cyllinder

% Circle geometry, selection of avoided sites
if strcmp(geometry, 'Circle')
    counter=1;
    for j=1:Ns
        if distance([position_matx(j,2), position_matx(j,3)])>=Radius;
            avoided_list(counter)=j;
            counter=counter+1;
        end
    end
    for state=1:Ns
        if sum(state==avoided_list)==1
            H(state,:)=0; H(:, state)=0;
        end
    end
end

% Plot of the finite geometry with names of the different position states
figure(1) 
plot(position_matx(:,2), position_matx(:,3), '.b', 'MarkerSize', 15)
for n=1:Ns                                                                   % Label states in the corresponding site
    if sum(n==avoided_list)==0
    name=num2str(n);
    text(position_matx(n,2), position_matx(n,3), name)
    end
end
for i=1:Ns
    for j=1:Ns
        
        if H(i,j)==1
            line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'r')
        else if H(i,j)==2
                line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'b')
            else if H(i,j)==3
                    line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'g')
                end
            end
        end
    end
end                                                                % Label corresponding couplings  (in order to work we need to set J1, J2, J3 different in H)

%% Driving process and eigenstates
% We now implement the three-step driving period

% Hamiltonians for each driving step
if strcmp( drive, 'Kitagawa')==1
    if strcmp( geometry, 'Strip')==1
         [HF, H1, H2, H3] = Floquet_kit(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
    else if strcmp( geometry, 'Cyllinder')==1
            H1=Finite_Hamiltonian_Cyllinder( delta, lambda*J_1, J_2, J_3, Nx, Ny, position_matx); % 1st step Kitagawa
            H2=Finite_Hamiltonian_Cyllinder( delta, J_1, lambda*J_2, J_3, Nx, Ny, position_matx); % 2nd step Kitagawa
            H3=Finite_Hamiltonian_Cyllinder( delta, J_1, J_2, lambda*J_3, Nx, Ny, position_matx); % 3rd step Kitagawa
            
        else if strcmp( geometry, 'Circle')==1
                H1=Finite_Hamiltonian( delta, lambda*J_1, J_2, J_3, Nx, Ny, position_matx); % 1st step Kitagawa
                H2=Finite_Hamiltonian( delta, J_1, lambda*J_2, J_3, Nx, Ny, position_matx); % 2nd step Kitagawa
                H3=Finite_Hamiltonian( delta, J_1, J_2, lambda*J_3, Nx, Ny, position_matx); % 3rd step Kitagawa
                
                count=0;
                for state=1:Ns
                    if sum(state==avoided_list)==1
                        H1(state-count,:)=[]; H1(:, state-count)=[];
                        H2(state-count,:)=[]; H2(:, state-count)=[];
                        H3(state-count,:)=[]; H3(:, state-count)=[];
                        H(state-count,:)=[]; H(:, state-count)=[];
                        position_matx(state-count,:)=[];
                        Psi_0(state-count)=[]; 
                        overlap_0(state-count)=[]; 
                        total_prob_0(state-count)=[];  
                        total_prob_pi(state-count)=[]; 
                        count=count+1;
                    end
                end
                Ns_real=length(H1(1,:));
                
            end
        end
    end
else if strcmp(drive, 'Rudner')==1
        if strcmp( geometry, 'Strip')==1
            H1=Finite_Hamiltonian( delta, J_1, 0, 0, Nx, Ny, position_matx); % 1st step Rudner
            H2=Finite_Hamiltonian( delta, 0, J_2, 0, Nx, Ny, position_matx); % 2nd step Rudner
            H3=Finite_Hamiltonian( delta, 0, 0, J_3, Nx, Ny, position_matx); % 3rd step Rudner
        else if strcmp( geometry, 'Cyllinder')==1
                H1=Finite_Hamiltonian_Cyllinder( delta, J_1, 0, 0, Nx, Ny, position_matx); % 1st step Rudner
                H2=Finite_Hamiltonian_Cyllinder( delta, 0, J_2, 0, Nx, Ny, position_matx); % 2nd step Rudner
                H3=Finite_Hamiltonian_Cyllinder( delta, 0, 0, J_3, Nx, Ny, position_matx); % 3rd step Rudner
            end
        end
    end
end

% %If we want to shift the origin of the strip
% position_matx(:,2)=position_matx(:,2)-Nx*3/8;
% position_matx(:,3)=position_matx(:,3)-y_sep*(Ny-1)/2;


% Impurities
if strcmp( impurities, 'On')==1
 H1(263,263)=50;
 H2(263,263)=50;
 H3(263,263)=50;
 H1(274,274)=50;
 H2(274,274)=50;
 H3(274,274)=50;
end

% % Stroboscopic time evolution (Floquet)
% U1=(expm(-1i*dt*H1))^n_t; U2=(expm(-1i*dt*H2))^n_t; U3=(expm(-1i*dt*H3))^n_t; % Evolution operators
% U=U3*U2*U1; HF=1i*logm(U);                                                    % Floquet Hamiltonian
[HF_states, HF_values]=eig(HF);                                               % Diagonal Floquet Hamiltonian
[eigenvalues, eigenstates]=order_eigenvalues(HF_values, HF_states);           % Ordered spectrum
quasienergy=diag(eigenvalues);                                                % Quasienergy of the finite strip

% Apply flattening
if strcmp(flattening, 'On')==1
quasienergy_flat=zeros(Ns,1);
quasienergy_flat(1:10)=linspace(-pi,-pi/2,10);
quasienergy_flat(11:Ns/2-8)=-pi/2;
quasienergy_flat(Ns/2-8:Ns/2+8)=linspace(-pi/2,pi/2,17);
quasienergy_flat(Ns/2+8:Ns-10)=pi/2;
quasienergy_flat(Ns-9:Ns)=linspace(pi/2,pi,10);
quasienergy=quasienergy_flat;
end

% Exact time evolution (through instantaneous hamiltonians)
[H1_states, H1_values]=eig(H1);                                               % Diagonal Hamiltonian
[eigenvalues1, eigenstates1]=order_eigenvalues(H1_values, H1_states);         % Ordered spectrum
energy1=diag(eigenvalues1);                                                   % Energy 
[H2_states, H2_values]=eig(H2);                                               % Diagonal Hamiltonian
[eigenvalues2, eigenstates2]=order_eigenvalues(H2_values, H2_states);         % Ordered spectrum
energy2=diag(eigenvalues2);                                                   % Energy 
[H3_states, H3_values]=eig(H3);                                               % Diagonal Hamiltonian
[eigenvalues3, eigenstates3]=order_eigenvalues(H3_values, H3_states);         % Ordered spectrum
energy3=diag(eigenvalues3);                                                   % Energy 



figure(2) % Plot the quasi energy spectrum
hold on
box on
plot(1:Ns_real, quasienergy, '.b', 'MarkerSize', 10)                          % Total quasienergy spectrum
% plot(1:12, quasienergy(1:12), '.g', 'MarkerSize', 10)                       % Pi-gap edge states  
% plot(Ns-11:Ns, quasienergy(Ns-11:Ns), '.g', 'MarkerSize', 10)               % 0-gap edge states  
% plot(Ns/2-12:Ns/2+12, quasienergy(Ns/2-12:Ns/2+12), '.r', 'MarkerSize', 10) % Pi-gap edge states  
%title('Quasienergy Spectrum')
xlabel('Eigenstate Number','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',25,'FontWeight','bold','Interpret','latex')
ylim([-pi, pi])
xlim([0 Ns_real])
ax = gca;
ax.FontSize =25; 
ax.XTick = [50,100,150,200,250,300];
%ax.XTickLabel = {'-\pi','-2','-1','0','1','2','\pi'};
ax.YTick = [-pi, -2, -1, 0, 1, 2, pi];
ax.YTickLabel = {'-\pi','-2','-1','0','1','2','\pi'};
print -depsc stripbands_[1,1].eps

%% Wave-packet Overlaps
% Next we initialise a wave-packet on the edge of the sample. We
% define the packet in the site basis, and then compute the evolution in
% the same basis using the Floquet eigenbasis. The packet may be gaussian
% using the gauss_weight function defined in the preamble, or could be made
% up of the certain states we want.

% Initial wave-packet
sign=1;
kick = 1;
for n=1:Ns_real
    
%%%%%%%%%%%%%%%%%%%%%%%% Gaussian packet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  if sum(n==avoided_list)==0;    
     if H(n,n)<0 %&& position_matx(n,3)>(y0-0.01) %(upper) (to alternate  same subl)   % Select a particular subalttice A>0 B<0
        position = [position_matx(n,2); position_matx(n,3)];                                        % Position of each state
        %kick = exp(-1i*K_BZ*position)*exp(1i*K_BZ*[x0; y0]);                                        % Momentum kick for the initial wavepcket
        Psi_0(n)=sign*kick*gauss_weight( position_matx(n,2),position_matx(n,3),sigma_x, sigma_y);   % Gaussian wavepacket in the sublattice
        sign = sign;                                                                                % Symmetric/ Antisymmetric wavepacket
%        else if H(n,n)>0 % && position_matx(n,3)>(y0-0.01) 
%           position = [position_matx(n,2); position_matx(n,3)];                                      % Position of each state
%           kick = exp(-1i*K2_BZ*position)*exp(1i*K_BZ*[x0; y0]); % Other sublattice
%           Psi_0(n)=sign*kick*gauss_weight( position_matx(n,2),position_matx(n,3),sigma_x, sigma_y); % Gaussian wavepacket in the sublattice
%            end
     end
  end

%%%%%%%%%%%%%%%%%%%%%%% Symmetry combinations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     if position_matx(n,3)>y0-0.01 %&& position_matx(n,2)>9.5 && position_matx(n,2)<23-9.5         % Select particular area of the lattice
%         if H(n,n)>0 % && position_matx(n,2)<12 %|| H(n,n)>0 && position_matx(n,2)>=18             % Select sublattice and mirror regions
%             position = [position_matx(n,2); position_matx(n,3)];                                  % Position of each state
%             % kick = exp(-1i*K_BZ*position)*exp(1i*K_BZ*[x0; y0]);  
%             Psi_0(n)=sign * kick;                                                                        % Support on state n
%             sign = -sign;                                                                         % Antisymmetric or symmetric wavefunction
% %         else if H(n,n)>0 && position_matx(n,2)>=12% && position_matx(n,2)<18                    % Other part of the mirror region               
% %                 Psi_0(n)=sign;                                                                  % Support on the state n                                     
% %            end
%         end
%     end
    
   
%%%%%%%%%%%%%%%%%%%%%%% Single site on the lattice %%%%%%%%%%%%%%%%%%%%%%%%
%    if n==2
%        Psi_0(n)=1;                                                                              % Support on state n
%    end

%     
%%%%%%%%%%%%%%%%%%%%%%% Superposition of edge states %%%%%%%%%%%%%%%%%%%%%%
%     if n>Ns-7 || n<=8
%             Psi_0=Psi_0+0.05*eigenstates(:,n);                                                  % Support on the set of eigenstates
%     
%     end
end
prob_density_0=abs(conj(Psi_0).*Psi_0);                                      % Initial probability density of the packet

% Overlap
for n=1:Ns_real 
    overlap_0(n)=abs(conj(transpose(eigenstates(:,n)))*Psi_0)^2;
end


figure(4)
% We use "object" programming to make a nice plot. We get three objects, hAx
% with the axis, hLine1 and 2 with the data points. We can modify them
% using object handlers as in Python adding .Handler. To set the properties
% of the axes is easier to use the set function.

%%%%%%%%%%%%%%%%%%%%%%%  MAIN FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
box on
hold on
[hAx,hLine1,hLine2] = plotyy(1:Ns_real, quasienergy, 1:Ns_real, conv_factor(pi, max(overlap_0))*overlap_0); 

set(gcf,'Position',[500 100 1300 1500])
pos = get(gca, 'Position');
set(gca, 'Position', [pos(1) pos(2)+0.05 pos(3) pos(4)-0.05]);

hLine1.Marker = '.';
hLine1.MarkerEdgeColor = 'b';
hLine1.MarkerSize = 10;
hLine2.LineStyle = '-';
hLine2.Color = 'r';

set(hAx(1),'ytick',[-pi, -pi/2, 0, pi/2, pi]);
set(hAx(1),'YTickLabel',{'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
set(hAx(2),'YTickLabel',[]);
set(hAx(1),'ylim',[-pi pi])
set(hAx(2),'ylim',[0 pi])
set(hAx, 'xlim', [0 Ns_real])
set(hAx(1),'ycolor','k')
set(hAx(2),'ycolor','k')

xlabel(hAx(1),'Eigenstate Number','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel(hAx(1),'$ET/J$','FontSize',25,'FontWeight','bold','Interpret','latex') % left y-axis
ylabel(hAx(2),'$\langle \psi_n \vert \psi_0 \rangle ^2 $ (arb. scale)','FontSize',25,'FontWeight','bold','Interpret','latex') % right y-axis

ax = gca;
ax.FontSize =25; 
%ax.XTick = [50,100,150,200,250,300];


%%%%%%%%%%%%%%%%%%%%%%%%%% INSET FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('Position',[0.15 0.8 0.5 0.07]) % We show the localisation in the overlap figure (left bottom width height)
hold on
box on

tol = 0.01;
aux_phase =[-pi:tol:pi];  % All possible phases a state may have
numPoints=length(aux_phase);
cmap = jet(numPoints);
for i=1:Ns_real
    
    if prob_density_0(i)<1e-15 % To avoid errors
        prob_density_0(i)=1e-15;
    end
    
    radius = conv_factor(50, max(prob_density_0))*prob_density_0(i);        % Radius proportional to the support on the sublattice
    phase = angle(Psi_0(i)); % Phase of the wavefunction on the particular site
    index = find(abs(aux_phase - phase)<=tol);                               % Index of the phase colormap
    color = cmap(index(1), :);                                              % Colour for the plot
    
    % plot(position_matx(i,2), position_matx(i,3), '.', 'Color', color, 'MarkerSize', radius)

    if real(Psi_0(i))>0 && position_matx(i,3)>y0-3*sqrt(3) % Plot probability densities (different color A and B)  (H(i, i)for sublattice plot)
        plot(position_matx(i,2), position_matx(i,3), '.m', 'MarkerSize', radius)
    else if real(Psi_0(i))<0 && position_matx(i,3)>y0-3*sqrt(3)
        plot(position_matx(i,2), position_matx(i,3), '.c', 'MarkerSize', radius)
        end 
    end
    
    % Honeycomb grid
    for j=1:Ns_real 
        
        if H(i,j)==1 %&& position_matx(i,3)>y0-3*sqrt(3)
            line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',1)
        else if H(i,j)==2 %&& position_matx(i,3)>y0-3*sqrt(3)
                line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
            else if H(i,j)==3 %&& position_matx(i,3)>y0-3*sqrt(3)
                    line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
                end
            end
        end
    end  
    
end
plot(x0, y0, '.k', 'MarkerSize', 8)  % Center of the wavepacket
xlim([-1 limx + 0.5])
ylim([limy-3 limy+0.5 ]) % upper edge
% ylim([-1 2]) % lower edge
%title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda), '    \sigma_x=', num2str(sigma_x), '    \sigma_y=', num2str(sigma_y)])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%set(gcf, 'PaperPositionMode', 'auto');
hold off;
xlabel('')
ylabel('')
%% Wave-packet evolution 
if strcmp( drive, 'Kitagawa')==1
    
   if strcmp( evolution, 'Exact')==1
            
            time=0:T/3:t_exact;
            for j=1:length(time)
                
                if j==1 % Initial step
                    Psi_t(:,j)=Psi_0;
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                end
                if step==1% Stepwise drive
                    Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy1, eigenstates1); % Evolution at time t, j many column vectors
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                    
                else if step==2
                        Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy2, eigenstates2); % Evolution at time t, j many column vectors
                        prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                        
                    else if step==3
                            Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy3, eigenstates3); % Evolution at time t, j many column vectors
                            prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                        end
                    end
                end
                
                step=step+1;
                if step>3
                    step=1;
                end % Change of step
            end
            
        else if strcmp( evolution, 'Stroboscopic')==1
                time=0:2*T:t_strob;                
                for j=1:length(time)
                    
                    Psi_t(:,j)= Wavepacket_evolution( Psi_0, time(j), quasienergy, eigenstates);   % Evolution at time t, j many column vectors
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));                           % Probability density of the packet
                    
                end
                
            end
        end
    
else if strcmp( drive, 'Rudner')==1
        
        if strcmp( evolution, 'Exact')==1
            
            time=0:T/3:t_exact;
            for j=1:length(time)
                
                if j==1 % Initial step
                    Psi_t(:,j)=Psi_0;
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                end
                if step==1% Stepwise drive
                    Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy1, eigenstates1); % Evolution at time t, j many column vectors
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                    
                else if step==2
                        Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy2, eigenstates2); % Evolution at time t, j many column vectors
                        prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                        
                    else if step==3
                            Psi_t(:,j)= Wavepacket_evolution( Psi_t(:, j-1), T/3, energy3, eigenstates3); % Evolution at time t, j many column vectors
                            prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                        end
                    end
                end
                
                step=step+1;
                if step>3
                    step=1;
                end % Change of step
            end
            
        else if strcmp( evolution, 'Stroboscopic')==1
                time=0:2*T:t_strob;                
                for j=1:length(time)
                    
                    Psi_t(:,j)= Wavepacket_evolution( Psi_0, time(j), quasienergy, eigenstates); % Evolution at time t, j many column vectors
                    prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
                    
                end
                
            end
        end
    end
end

% Movie of the evolution
if strcmp(movie, 'On')

figure(3) % Movie of the evolving wavepacket
% Initialise movie
framePerSec=2;
filenamefig=[pwd '/wavepacket_evolution1'];
vSim=VideoWriter([filenamefig '.wav']); % Creates a video archive in the current folder
vSim.FrameRate=framePerSec;
open(vSim);

for s=1:length(time) % Movie of the evolution
    
   
    prob_density_0(prob_density_0<1e-15)=1e-15;  % To avoid errors
    
    
    if strcmp( evolution, 'Stroboscopic')==1% Title selection
        heading=['\omega/J=', num2str(omega), '    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda), ...
            '    \sigma_x=', num2str(sigma_x), '    \sigma_y=', num2str(sigma_y),'    t=', num2str(s-1),'T'];
    else if strcmp( evolution, 'Exact')==1
            heading=['\omega/J=', num2str(omega), '    \Delta/J=', num2str(delta), ...
                '    \sigma_x=', num2str(sigma_x), '    \sigma_y=', num2str(sigma_y),'    t=', num2str(s-1),'T/3'];
        end
    end 
    
    hold on
    box on
      
    for i=1:Ns_real    % Plot evolution
        
        %if sum(i==avoided_list)==0
        
        if H(i, i)>0 % Plot probability densities (different color A and B)  
            plot(position_matx(i,2), position_matx(i,3), '.r', 'MarkerSize', conv_factor(100, max(prob_density(:,s)))*prob_density(i,s))
        else           
            plot(position_matx(i,2), position_matx(i,3), '.b', 'MarkerSize', conv_factor(100, max(prob_density(:,s)))*prob_density(i,s))
        end   
        if H1(i,i)>5 % Signal Impurities
            plot(position_matx(i,2), position_matx(i,3), '.k', 'MarkerSize', 30)
        end
        for j=1:Ns_real % Honeycomb grid    
            
            if H(i,j)==1
                line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',1)
            else if H(i,j)==2
                    line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
                else if H(i,j)==3
                        line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
                    end
                end
            end
        end  
        %end
        
    end 

    xlim([-1 limx])
    ylim([-1 limy])
    title(heading);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    set(gcf,'Position',[900 100 1000 1000])  
    
    frame = getframe(gcf); % Gets the the plot as the frame for the video 
    writeVideo(vSim,frame);
    hold off;
    
    
    clf % Clear after each iteration
    
end 
close(vSim); % Close the animation


end


% print -depsc try2.eps

%% Total probability of each gap
% % Here we add all the probabilities from the different states belonging to
% % the same gap. It is IMPORTANT TO SELECT THE NUMBER OF EDGE STATES WE TAKE
% % AS THEY VARY FROM PHASE TO PHASE AND FOR DIFFERENT GEOMETRIES
% 
% % Zero gap
% total_prob_0=0;
% for n=Ns/2-7:Ns/2+7
% total_prob_0=total_prob_0+abs(conj(eigenstates(:,n)).*eigenstates(:,n));
% end
% % Pi gap
% total_prob_pi=0;
% for n=[1:8,Ns-8:Ns] 
% total_prob_pi=total_prob_pi+abs(conj(eigenstates(:,n)).*eigenstates(:,n));
% end
% 
% figure(6)
% box on
% hold on
% for i=1:Ns
%     
%     if H(i,i)>0 % Plot probability densities (different color A and B)
%         plot(position_matx(i,2), position_matx(i,3), '.r', 'MarkerSize', conv_factor(100, max(total_prob_0))*total_prob_0(i))
%     else
%         plot(position_matx(i,2), position_matx(i,3), '.b', 'MarkerSize', conv_factor(100, max(total_prob_0))*total_prob_0(i))
%     end 
%     for j=1:Ns
%         
%         if H(i,j)==1
%             line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',1)
%         else if H(i,j)==2
%                 line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%             else if H(i,j)==3
%                     line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%                 end
%             end
%         end
%     end  % Honeycomb grid
%     
% end
% xlim([-1 limx])
% ylim([-1 limy])
% %title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda)])
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gcf,'Position',[700 400 1000 600])
% hold off;
% print -depsc position0_[1,0].eps
% 
% figure(7)
% box on
% hold on
% for i=1:Ns
%     
%     if H(i,i)>0 % Plot probability densities (different color A and B)
%         plot(position_matx(i,2), position_matx(i,3), '.r', 'MarkerSize', conv_factor(100, max(total_prob_pi))*total_prob_pi(i))
%     else
%         plot(position_matx(i,2), position_matx(i,3), '.b', 'MarkerSize', conv_factor(100, max(total_prob_pi))*total_prob_pi(i))
%     end 
%     for j=1:Ns % Honeycomb grid
%         
%         if H(i,j)==1
%             line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',1)
%         else if H(i,j)==2
%                 line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%             else if H(i,j)==3
%                     line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%                 end
%             end
%         end
%     end  
%     
% end
% xlim([-1 limx])
% ylim([-1 limy])
% %title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda)])
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gcf,'Position',[700 400 1000 600])
% hold off;
% %print -depsc positionpi_[0,1].eps

%% Symmetry analysis
% Sign of the initial state selection. This only make sense when we enter
% the selection by hand and it's real.

% figure(8) % Snapshot of t=0 in the evolution
% hold on
% box on
% for i=1:Ns_real
%     
%     if prob_density_0(i)<1e-15 % To avoid errors
%         prob_density_0(i)=1e-15;
%     end
%     
%     if Psi_0(i)>0 % Plot probability densities with sign color of the amplitude 
%         plot(position_matx(i,2), position_matx(i,3), '.m', 'MarkerSize', conv_factor(100, max(prob_density_0))*prob_density_0(i))
%     else if Psi_0(i)<0
%         plot(position_matx(i,2), position_matx(i,3), '.c', 'MarkerSize', conv_factor(100, max(prob_density_0))*prob_density_0(i))
%     end 
%     end
%     for j=1:Ns_real % Honeycomb grid
%         
%         if H(i,j)==1
%             line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',1)
%         else if H(i,j)==2
%                 line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%             else if H(i,j)==3
%                     line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',1)
%                 end
%             end
%         end
%     end  
%     
% end
% xlim([-1 limx])
% ylim([-1 limy])
% title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda), '    \sigma_x=', num2str(sigma_x), '    \sigma_y=', num2str(sigma_y)])
% set(gca,'YTickLabel',[]);
% set(gca,'XTickLabel',[]);
% set(gcf,'Position',[900 400 1000 600])
% hold off;

% %% Evolution figures
% % We make figures to show the evolution of the wave packets without needing
% % to play a movie. Inteded for publishing results.
% 
% cont=1;
% for n=1:Ns_real % We take the probability density of sites at the edge,
%            % and take it into a matrix as they evolve in time
%     
% if position_matx(n,3)>y0-sqrt(3)
%     
%     evolution_matx(:,cont)=prob_density(n,:);
%     cont=cont+1;
% end
%    
% 
% end
% 
% close all
% figure(9)
% s=surf(evolution_matx);
% set(gcf,'Position',[900 400 600 500])
% box on
% view(2)
% c=jet;
% colormap(c(10:50,:));
% caxis([0 2e-3])
% colorbar
% xlabel('Edge Site','FontSize',25,'FontWeight','bold','Interpret','latex')
% ylabel('Time [T]','FontSize',25,'FontWeight','bold','Interpret','latex')
% h = colorbar;
% set(get(h,'label'),'string','$\vert \Psi \vert^2$ [arbitrary units]','FontSize',20,'FontWeight','bold','Interpret','latex');
% xlim([1 32])
% ylim([1 length(time)])
% s.EdgeColor = 'flat'
% ax = gca;
% ax.FontSize =25; 
% ax.XTick = [1,10,20,30];
% 
% print -depsc uppercurrent_[-1,1].eps
% 
% %% Graphene sheet
% % This is just to produce figures of thehoneycomb grid we are using
% close all
% figure(10)
% hold on
% box on
% 
% for i=1:Ns
%       
%     for j=1:Ns
%         
%         if H(i,j)==1
%             line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k','LineWidth',2)
%         else if H(i,j)==2
%                 line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',2)
%             else if H(i,j)==3
%                     line([position_matx(i,2) position_matx(j,2)], [position_matx(i,3) position_matx(j,3)], 'Color', 'k', 'LineWidth',2)
%                 end
%             end
%         end
%     end  % Honeycomb grid
%     
% end
% for i=1:Ns
% if H(i,i)>0 % Plot probability densities (different color A and B
%         
%         scatter(position_matx(i,2), position_matx(i,3),100, 'b', 'filled','MarkerEdgeColor','k')
%     else
%         scatter(position_matx(i,2), position_matx(i,3),100, 'r', 'filled','MarkerEdgeColor','k')
% end 
% end
% %  text(0.9,2,'$J_1$','FontSize',17,'FontWeight','bold','Interpret','latex')
% %  text(0.9,-0.3,'$ J_1$','FontSize',17,'FontWeight','bold','Interpret','latex')
% % text(1.8,1.4,'$J_2$','FontSize',17,'FontWeight','bold','Interpret','latex')
% % text(-0.1,0.3,'$ J_2$','FontSize',17,'FontWeight','bold','Interpret','latex')
% % text(1.8,0.3,'$\lambda J_3$','FontSize',17,'FontWeight','bold','Interpret','latex')
% % text(-0.2,1.35,'$\lambda J_3$','FontSize',17,'FontWeight','bold','Interpret','latex')
% % text(0.8,0.85,'$(3.)$','FontSize',20,'FontWeight','bold','Interpret','latex')
% text(2,-0.5,'$J_1$','FontSize',17,'FontWeight','bold','Interpret','latex')
%  text(1.5,-0.5,'$J_2$','FontSize',17,'FontWeight','bold','Interpret','latex')
%  text(0.5,-0.5,'$J_3$','FontSize',17,'FontWeight','bold','Interpret','latex')
%  text(0.1,-0.5,'$A \quad B$','FontSize',17,'FontWeight','bold','Interpret','latex')
%   text(2.5,-0.5,'$\overrightarrow{a}_1$','FontSize',17,'FontWeight','bold','Interpret','latex')
%     text(4,-0.5,'$\overrightarrow{a}_2$','FontSize',17,'FontWeight','bold','Interpret','latex')
% 
% xlim([-1 limx])
% ylim([-1 limy])
% set(gcf,'Position',[900 400 600 400])
% 
%% Symmetry connection between different gaps
% 
% for index=1:Ns % Taking the (centered) branchcut shift for the strip  Floquet Hamiltonian
%    
%     if eigenvalues(index, index)>0
%       HF_pi = HF_pi + real(eigenvalues(index,index))* eigenstates(:,index)*(eigenstates(:,index)');
%       HF_0 = HF_0 + real(eigenvalues(index,index)-pi)* eigenstates(:,index)*(eigenstates(:,index)');
%     else
%       HF_pi = HF_pi + real(eigenvalues(index,index))* eigenstates(:,index)*(eigenstates(:,index)');
%       HF_0 = HF_0 + real(eigenvalues(index,index)+pi)* eigenstates(:,index)*(eigenstates(:,index)');
%     end  
% end
% 
% [aux1, aux2]=eig(HF_0);
% [aux3, aux4]=order_eigenvalues(aux2, aux1); % Ordered spectrum
% aux5=diag(aux3); % Quasienergy of the finite strip
% 
% figure(11)
% plot(1:Ns, aux5, '.b')
% xlim([0 Ns_real])
% 
% v1=ones(1, Ns/2)*(-1i);
% v2=ones(1, Ns/2)*(1i);
% v3=diag(v1,Ns/2);
% v4=diag(v2,-Ns/2);
% 
% sigmay_strip=zeros(Ns)+v3+v4;
% 
% symmetry1 = sigmay_strip * HF_pi * sigmay_strip;
% symmetry2 = conj(HF_0);
% 

