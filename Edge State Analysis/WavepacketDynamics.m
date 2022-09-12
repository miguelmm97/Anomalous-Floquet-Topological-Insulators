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

%% Variables

% Code specification
drive='Rudner';                                                           % Select driving approach
geometry= 'Cyllinder';% 'Cyllinder' % 'Circle'                              % Select geometry
evolution= 'Exact'; % 'Exact'                                        % Evolution for the Rudner drive
flattening= 'Off'; % Off                                                    % Apply flattening to the quasienergy spectrum
impurities= 'Off';                                                          % Put impurities in the sample
movie= 'Off';                                                               % If we want to produce a movie of the evolution

% Region in the phase diagram
omega=4/3;     % 6.5 7 5 (Sweet spot Rudner drive w=4/3)                    % Frequency of the drive ([J])
delta_plot=0.5;
delta=0; % 2 3                                                              % Potential offset ([J])
lambda=3;   % 2.5 4.5 3.4                                                   % Degree of anisotropy
J=1;                                                                        % Coupling constant
J_1=J; J_2=J; J_3=J;                                                        % Hoping amplitudes ([J])
bc = 0;                                                                     % branchcut at 0 if bc=1

% Time
T=2*pi/omega;                                                               % Period
t0 = 0*T/20;                                                               % Initial time
t_strob=60*T; t_exact=60*T;                                                 % Final time of the packet evolution

% Geometry
Nx=104; % 32 Necessarilly mod(Nx,4)=0 for zig-zag edge along y              % Number of x sites
Ny=25; % 15 Necessarilly odd for armchair edge along x                      % Number of y sites
Ns=Nx*Ny/2;   Ns_real=Ns; %(Ns_real only changes for circular geometry)     % Total number of states
y_sep=sqrt(3)/2;    x_sep=3;                                                % Separation between y layers 
layers_x=Nx*3/8;                                                            % x length of each cell
center=[Nx*3/8, y_sep*(Ny-1)/2];                                            % Center of the circle geometry
Radius=Nx*3/8;                                                              % Radius of the circle geometry

% Brillouin Zone
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);
gamma_BZ=[0,0];
K_BZ=[l,0];
K2_BZ=[f,s];
nodal_line=1.01;

% Wavepackets
x01=((Nx/4)*x_sep/2)-1.5; y01=0; %y_sep*(Ny-1) ;    %0                       % Initial position of the wave-packet's center 1
x02=((Nx/4)*x_sep/2)-1.5; y02=0; %y_sep*(Ny-1) ;     %0                      % Initial position of the wave-packet's center 2
sigma_x1=1; sigma_y1=0.5;                                                   % Widths 1
sigma_x2=1; sigma_y2=0.5;                                                   % Widths 2
symmetry=1;                                                                 % Symmetric/ Antisymmetric wavepacket
sign=1;                                                                     % Sign of the wavepacket
mom_kick =[0,0]; %Z; %K_BZ/2;                                               % Momentum kick
edge = 'Lower';

% Declarations
position_matx=zeros(Ns,3);                                                   % 1st colum state number, 2nd x position, 3rd y position
Psi_0=zeros(Ns,1);                                                           % Initial wavepacket
overlap_0=zeros(Ns,1);                                                       % Initial overlap of the wave packet
total_prob_0=zeros(Ns,1);                                                    % Total prob for 0 gap
total_prob_pi=zeros(Ns,1);                                                   % Total prob pi gap
step=0;                                                                      % Step counter in the evolution of the Rudner drive
avoided_list=0;                                                              % Thrown away points in the circle geometry

% Functions
gauss_weight=@(x, y, x0, y0, sigma_x, sigma_y) ... 
           exp(-((x-x0)^2)/(4*sigma_x^2))* exp(-((y-y0)^2)/(4*sigma_y^2));   % Gaussian probability amplitude (not normalised)
conv_factor=@(value, maximum) value/maximum;                                 % Conversion factor to plot probability densities
distance=@(vector) sqrt((vector(1)-center(1))^2+(vector(2)-center(2))^2);    % Norm function

%% State-position assignment and Hamiltonian Construction
disp('Calculating geometry...')

% Construction of the position matrix (n_state, xpos, ypos)
position_matx=position_matrix(Nx, Ny);                                       % State position and labelling
limx=max(position_matx(:,2))+1;  limy=max(position_matx(:,3))+1;             % Usual limits of the strip to plot later on
H=Finite_Hamiltonian( delta_plot, 1, 2, 3, Nx, Ny, position_matx);           % Test hamiltonian strip

% Circle geometry, selection of avoided sites
switch geometry
    case 'Circle'
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

%% Micromotion and stroboscopic Hamiltonians and spectrum
disp('Calculating Hamiltonians...')

% Hamiltonians for each driving step
switch drive   
   
    case 'Kitagawa'          
        switch geometry
            case 'Strip'
                [HF, H1, H2, H3] = Floquet_kit_strip(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, T, t0);
            case 'Cyllinder'
                [HF, H1, H2, H3] = Floquet_kit_cyl(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, T, t0);
            case 'Circle'
                 [HF, H1, H2, H3] = Floquet_kit_strip(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);  
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
            otherwise
                 disp('Unknown geometry!')
        end
        
    % Rudner drive
   
    case 'Rudner'       
       switch geometry
           case 'Strip'
               [HF, H1, H2, H3] = Floquet_rud_strip(delta, J_1, J_2, J_3, Nx, Ny, position_matx, T, t0);
           case 'Cyllinder'
                [HF, H1, H2, H3] = Floquet_rud_cyl(delta, J_1, J_2, J_3, Nx, Ny, position_matx, T, t0);
           case 'Circle'
                 [HF, H1, H2, H3] = Floquet_rud_strip(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
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
           otherwise
                 disp('Unknown geometry!')
       end
       
    % Raise   
    otherwise
        disp('Unknown drive!')
end

% Impurities
switch impurities
    case 'On'
      H1(263,263)=50;
      H2(263,263)=50;
      H3(263,263)=50;
      H1(274,274)=50;
      H2(274,274)=50;
      H3(274,274)=50;
    otherwise
end

% Stroboscopic and micromotion spectrum
% qen = diag(qen)
% if max(abs(imag(quasienergy(:)))) > 1e-6, warning("Imaginary qen")
% [quasienergy, eigenstates] = eig(HF);
% [quasienergy, ind] = sort(real(qen));
% vec = vec(:,ind);

[quasienergy, eigenstates]=spectrum(HF);                                      % Spectrum Floquet Hamiltonian
[energy1, eigenstates1]=spectrum(H1);                                         % Spectrum 1st Hamiltonian
[energy2, eigenstates2]=spectrum(H2);                                         % Spectrum 2nd Hamiltonian
[energy3, eigenstates3]=spectrum(H3);                                         % Spectrum 3rd Hamiltonian
quasienergy(quasienergy<0) = quasienergy(quasienergy<0)+ 2*pi*bc;             % Specified branchcut

% Apply flattening
switch flattening
    case 'On'
        quasienergy_flat=zeros(Ns,1);
        quasienergy_flat(1:10)=linspace(-pi,-pi/2,10);
        quasienergy_flat(11:Ns/2-8)=-pi/2;
        quasienergy_flat(Ns/2-8:Ns/2+8)=linspace(-pi/2,pi/2,17);
        quasienergy_flat(Ns/2+8:Ns-10)=pi/2;
        quasienergy_flat(Ns-9:Ns)=linspace(pi/2,pi,10);
        quasienergy=quasienergy_flat;
end
      
%% Initial Wave-packet and Overlaps
disp('Calculating initial wavepacket...')

% Initial wave-packet
switch geometry
    
    case {'Strip', 'Cyllinder'}        
        for n=1:Ns_real 
            % A sublattice
             if H(n,n)>0 % && position_matx(n,3)<y01+0.1 ; % y01-sqrt(3)/2+0.1
                position = [position_matx(n,2); position_matx(n,3)];                                        % Position of each state
                kick = exp(-1i*mom_kick*position); %*exp(1i*mom_kick*[x0; y0]);                             % Momentum kick for the initial wavepcket
                Psi_0(n)=sign*kick*(gauss_weight(position_matx(n,2),position_matx(n,3), x01, y01, sigma_x1, sigma_y1)+... 
                gauss_weight(position_matx(n,2), position_matx(n,3), x02, y02, sigma_x2, sigma_y2));        % Gaussian wavepacket in the sublattice               
                sign = symmetry*sign;                                                                       % Symmetric/ Antisymmetric wavepacket       
          
            % B sublattice
            else if H(n,n)<0 % && position_matx(n,3)<y01 +0.1; % y01+sqrt(3)/2 -0.1
                    position = [position_matx(n,2); position_matx(n,3)];                                        % Position of each state
                    kick = exp(-1i*mom_kick*position); %*exp(1i*mom_kick*[x0; y0]);                             % Momentum kick for the initial wavepcket
                    Psi_0(n)=sign*kick*(gauss_weight( position_matx(n,2),position_matx(n,3), x01, y01, sigma_x1, sigma_y1)+... 
                    gauss_weight( position_matx(n,2),position_matx(n,3), x02, y02, sigma_x2, sigma_y2));        % Gaussian wavepacket in the sublattice     
                    sign = symmetry*sign;                                                                       % Symmetric/ Antisymmetric wavepacket
                end
            end
        end
        
    case 'Circle'
        disp('Must specify initial wavepakcet!')
        
    otherwise
        disp('Wrong geometry!')
end 

% Overlap and initial distribution
% Psi_0 = eigenstates(:,390);                                              % Select a particular eigenstate
norm_psi0 = sqrt(transpose(conj(Psi_0))*Psi_0);                            % Normalisation
Psi_0 = Psi_0 / norm_psi0;                                                 % Normalisation
prob_density_0=abs(conj(Psi_0).*Psi_0);                                    % Initial probability density of the packet
for n=1:Ns_real                                                            % Initial overlap
    overlap_0(n)=abs(conj(transpose(eigenstates(:,n)))*Psi_0)^2;
end

%% Wave-packet evolution 
disp('Calculating evolution...')  

switch evolution
    case 'Exact'  
        step=0;
        time=0:T/3:t_exact;
        Psi_t=zeros(Ns,length(time));
        for j=1:length(time)
            % Initial step
            if j==1
                Psi_t(:,j)=Psi_0;
                prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
            end
            % Stepwise drive
            if step==1
                coefs1 = eigenstates1'*Psi_t(:, j-1);
                coefs_t1 = coefs1.*exp(-1i*energy1*T/3);                  % Initial overlap with eigenstates
                Psi_t(:,j) = eigenstates1*coefs_t1;                               % Wavefunction at t
                prob_density(:,j)=abs(Psi_t(:,j)).^2;                           % Probability at t
            else if step==2
                    coefs2 = eigenstates2'*Psi_t(:, j-1);
                    coefs_t2 = coefs2.*exp(-1i*energy2*T/3);                  % Initial overlap with eigenstates
                    Psi_t(:,j) = eigenstates2*coefs_t2;                               % Wavefunction at t
                    prob_density(:,j)=abs(Psi_t(:,j)).^2;     
                else if step==3
                        coefs3 = eigenstates3'*Psi_t(:, j-1);
                        coefs_t3 = coefs3.*exp(-1i*energy3*T/3);                  % Initial overlap with eigenstates
                        Psi_t(:,j) = eigenstates3*coefs_t3;                               % Wavefunction at t
                        prob_density(:,j)=abs(Psi_t(:,j)).^2;     
                    end
                end
            end
            
            % Change of step
            step=step+1;
            if step>3
                step=1;
            end
        end
        
    case 'Stroboscopic'
        time=0:T:t_strob;
        Psi_t=zeros(Ns,length(time));
        coefs = eigenstates'*Psi_0;
        for j=1:length(time)   
            coefs_t = coefs.*exp(-1i*quasienergy*time(j));                  % Initial overlap with eigenstates
            Psi_t(:,j) = eigenstates*coefs_t;                               % Wavefunction at t
            prob_density(:,j)=abs(Psi_t(:,j)).^2;                           % Probability at t
        end
        
    otherwise
        disp('Unknown evolution!')
end

% We take the probability density of sites at the edge
delta_x = 1/2;                                                                    % Size of the batch for x disctrtisation
evolution_matx = Edge_evolution(prob_density, position_matx, Ny, delta_x, edge);  % Evolution only at the edge

%% Saving data for the averaged dynamics
% str1 = 't0'; str2 = num2str(100*round(t0,2)); str3 = '.mat';
% file = strcat(str1, str2, str3);
% cd Averaged_dynamics
% switch edge
%     case 'Upper'        
%        cd Upper_edge_data % upper edge
%     case 'Lower'
%        cd Lower_edge_data % lower edge
% end
% save(file, 'evolution_matx')
% return

%% Figures

% GEOMETRY WITH LABELLING
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
title('Geometry and site labeling','FontSize',25,'FontWeight','bold','Interpret','latex')

% QUASIENERGY BANDS
figure(2) 
hold on
box on
plot(1:Ns_real, quasienergy.*T, '.b', 'MarkerSize', 10)                          % Total quasienergy spectrum
% plot(1:12, quasienergy(1:12), '.g', 'MarkerSize', 10)                       % Pi-gap edge states  
% plot(Ns-11:Ns, quasienergy(Ns-11:Ns), '.g', 'MarkerSize', 10)               % 0-gap edge states  
% plot(Ns/2-12:Ns/2+12, quasienergy(Ns/2-12:Ns/2+12), '.r', 'MarkerSize', 10) % Pi-gap edge states  
%title('Quasienergy Spectrum')
xlabel('Eigenstate Number','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',25,'FontWeight','bold','Interpret','latex')
title([drive, ' $\omega= $', num2str(omega), ', $\delta= $', num2str(delta), ', $\lambda= $', num2str(lambda),...
     ', $N_x=$ ', num2str(Nx), ', $N_y=$ ', num2str(Ny),', $t_0=$ ', num2str(t0)],'FontSize',20,'FontWeight','bold','Interpret','latex'); 
switch bc
    case 1
ylim([-2*pi, 2*pi])
    case 0
ylim([-pi, pi])   
end 

xlim([0 Ns_real])
ax = gca;
ax.FontSize =25; 
% ax.XTick = [50,100,150,200,250,300];
%ax.XTickLabel = {'-\pi','-2','-1','0','1','2','\pi'};
ax.YTick = [-pi, -2, -1, 0, 1, 2, pi];
ax.YTickLabel = {'-\pi','-2','-1','0','1','2','\pi'};
% print -depsc stripbands_[1,1].eps



% INITIAL OVERLAP
figure(3)
% Main figure
box on
hold on
[hAx,hLine1,hLine2] = plotyy(1:Ns_real, quasienergy.*T, 1:Ns_real, conv_factor(pi, max(overlap_0))*overlap_0); 

set(gcf,'Position',[400 20 1200 850])
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
switch bc
    case 1
set(hAx(1),'ylim',[-2*pi 2*pi])
    case 0
set(hAx(1),'ylim',[-pi pi])  
end 
set(hAx(2),'ylim',[0 pi])
set(hAx, 'xlim', [0 Ns_real])
set(hAx(1),'ycolor','k')
set(hAx(2),'ycolor','k')

xlabel(hAx(1),'Eigenstate Number','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel(hAx(1),'$ET/J$','FontSize',25,'FontWeight','bold','Interpret','latex') % left y-axis
ylabel(hAx(2),'$\langle \psi_n \vert \psi_0 \rangle ^2 $ (arb. scale)','FontSize',25,'FontWeight','bold','Interpret','latex') % right y-axis
title([drive, ' $\omega= $', num2str(omega), ', $\delta= $', num2str(delta), ', $\lambda= $', num2str(lambda),...
     ', $N_x=$ ', num2str(Nx), ', $N_y=$ ', num2str(Ny),', $t_0=$ ', num2str(t0)],'FontSize',20,'FontWeight','bold','Interpret','latex'); 
ax = gca;
ax.FontSize =25; 
%ax.XTick = [50,100,150,200,250,300];


% Inset figure
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

    if real(Psi_0(i))>0 && position_matx(i,3)>y01-3*sqrt(3) % Plot probability densities (different color A and B)  (H(i, i)for sublattice plot)
        plot(position_matx(i,2), position_matx(i,3), '.m', 'MarkerSize', radius)
    else if real(Psi_0(i))<0 && position_matx(i,3)>y01-3*sqrt(3)
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
plot(x01, y01, '.k', 'MarkerSize', 16)  % Center of the wavepacket
plot(x02, y02, '.k', 'MarkerSize', 16)  % Center of the wavepacket
xlim([-1 limx + 0.5])
switch edge
    case 'Upper'        
       ylim([limy-4 limy+0.5 ]) % upper edge
    case 'Lower'
       ylim([-1 2]) % lower edge
end
%title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda), '    \sigma_x=', num2str(sigma_x), '    \sigma_y=', num2str(sigma_y)])
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
%set(gcf, 'PaperPositionMode', 'auto');
hold off;
xlabel('')
ylabel('')


% EVOLUTION FIGURE
figure(4)
s=surf(evolution_matx);
set(gcf,'Position',[700 200 1000 700])
box on
view(2)
c=jet;
colormap(c);
caxis([0 0.06])
colorbar
xlabel('$x[a]$','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel('Time [T]','FontSize',25,'FontWeight','bold','Interpret','latex')
title([drive, ' $\omega= $', num2str(omega), ', $\delta= $', num2str(delta), ', $\lambda= $', num2str(lambda),...
     ', $N_x=$ ', num2str(Nx), ', $N_y=$ ', num2str(Ny),', $t_0=$ ', num2str(t0)],'FontSize',20,'FontWeight','bold','Interpret','latex'); 
h = colorbar;
set(get(h,'label'),'string','$\vert \Psi \vert^2$ [arbitrary units]','FontSize',20,'FontWeight','bold','Interpret','latex');
xlim([1 length(evolution_matx(1, :))])
ylim([1 length(time)])
s.EdgeColor = 'flat'
ax = gca;
ax.FontSize =25; 
ax.XTick = [0:10:length(evolution_matx(1, :))];




% MOVIE OF THE EVOLUTION
if strcmp(movie, 'On')

figure(5) % Movie of the evolving wavepacket
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
            '    \sigma_x=', num2str(sigma_x1), '    \sigma_y=', num2str(sigma_y1),'    t=', num2str(s-1),'T'];
    else if strcmp( evolution, 'Exact')==1
            heading=['\omega/J=', num2str(omega), '    \Delta/J=', num2str(delta), ...
                '    \sigma_x=', num2str(sigma_x1), '    \sigma_y=', num2str(sigma_y1),'    t=', num2str(s-1),'T/3'];
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
    set(gcf,'Position',[500 50 1100 700])  
    
    frame = getframe(gcf); % Gets the the plot as the frame for the video 
    writeVideo(vSim,frame);
    hold off;
    
    
    clf % Clear after each iteration
    
end 
close(vSim); % Close the animation


end






