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
drive='Rudner';                                                           % 'Kitagawa' % Select driving approach
geometry= 'Strip';% 'Cyllinder' % 'Circle'                              % Select geometry
evolution= 'Stroboscopic'; % 'Exact'                                        % Evolution for the Rudner drive
flattening= 'Off'; % Off                                                    % Apply flattening to the quasienergy spectrum
impurities= 'Off';                                                          % Put impurities in the sample
movie= 'Off';                                                               % If we want to produce a movie of the evolution

% Region in the phase diagram
omega=5;     % 6.5 7 5 (Sweet spot Rudner drive w=4/3)                      % Frequency of the drive ([J])
delta_plot=0.5;
delta=2; % 2 3                                                              % Potential offset ([J])
lambda=3.4;   % 2.5 4.5 3.4                                                 % Degree of anisotropy
J=1;                                                                        % Coupling constant
J_1=J; J_2=J; J_3=J;                                                        % Hoping amplitudes ([J])

% Geometry
Nx=32; % Necessarilly mod(Nx,4)=0 for zig-zag edge along y                  % Number of x sites
Ny=15; % Necessarilly odd for armchair edge along x                         % Number of y sites
Ns=Nx*Ny/2;   Ns_real=Ns;                                                   % Total number of states
y_sep=sqrt(3)/2;                                                            % Separation between y layers so that every site is equidistant from the others in units of a=1
x_sep=3;  layers_x=Nx*3/8;                                                  % x length of each cell
center=[Nx*3/8, y_sep*(Ny-1)/2];                                            % Center of the circle geometry
Radius=Nx*3/8;                                                              % Radius of the circle geometry

% Wavepackets
x0=((Nx/4)*x_sep/2)-0.5; y0=y_sep*(Ny-1);                                % Initial position of the wave-packet's center
sigma_x=2; sigma_y=1;                                                       % Widths
edge = 'Upper';

%% Definitions

% Time
n_t=100;                                                                    % Time discretisation
T=2*pi/omega;                                                               % Period
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];      % Time grid
dt=t(2)-t(1);                                                               % Time step
t_strob=60*T; t_exact=30*T;                                                 % Final time of the packet evolution

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

% Declarations
position_matx=zeros(Ns,3);                                                   % 1st colum state number, 2nd x position, 3rd y position
Psi_0=zeros(Ns,1);                                                           % Initial wavepacket
overlap_0=zeros(Ns,1);                                                       % Initial overlap of the wave packet
total_prob_0=zeros(Ns,1);                                                    % Total prob for 0 gap
total_prob_pi=zeros(Ns,1);                                                   % Total prob pi gap
step=0;                                                                      % Step counter in the evolution of the Rudner drive
avoided_list=0;                                                              % Thrown away points in the circle geometry
HF_pi = zeros(Ns); HF_0=zeros(Ns);                                           % Definition for hamiltonians at both gaps for symmetry connection


% Functions
gauss_weight=@(x, y, sigma_x, sigma_y) ... 
           exp(-((x-x0)^2)/(4*sigma_x^2))* exp(-((y-y0)^2)/(4*sigma_y^2));  
conv_factor=@(value, maximum) value/maximum;                                 % Conversion factor to plot probability densities
distance=@(vector) sqrt((vector(1)-center(1))^2+(vector(2)-center(2))^2);    % Norm function


%% State-position assignment and Hamiltonian Construction

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

% Hamiltonians for each driving step
switch drive   
   
    case 'Kitagawa'          
        switch geometry
            case 'Strip'
                [HF, H1, H2, H3] = Floquet_kit_strip(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
            case 'Cyllinder'
                [HF, H1, H2, H3] = Floquet_kit_cyl(delta, lambda, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
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
               [HF, H1, H2, H3] = Floquet_rud_strip(delta, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
           case 'Cyllinder'
                [HF, H1, H2, H3] = Floquet_rud_cyl(delta, J_1, J_2, J_3, Nx, Ny, position_matx, dt, n_t);
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
[quasienergy, eigenstates]=spectrum(HF);                                      % Spectrum Floquet Hamiltonian
[energy1, eigenstates1]=spectrum(H1);                                         % Spectrum 1st Hamiltonian
[energy2, eigenstates2]=spectrum(H2);                                         % Spectrum 2nd Hamiltonian
[energy3, eigenstates3]=spectrum(H3);                                         % Spectrum 3rd Hamiltonian

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

      
%% Wave-packet Overlaps

% Initial wave-packet
symmetry=1;
sign=1;
mom_kick =[0,0];%nodal_line*[cos(pi/2), sin(pi/2)];
switch geometry
    
    case {'Strip', 'Cyllinder'}        
        for n=1:Ns_real 
            % A sublattice
            if H(n,n)>0 %&& position_matx(n,3)>y0-sqrt(3)
                position = [position_matx(n,2); position_matx(n,3)];                                        % Position of each state
                kick = exp(-1i*mom_kick*position)*exp(1i*mom_kick*[x0; y0]);                                % Momentum kick for the initial wavepcket
                Psi_0(n)=sign*kick*gauss_weight( position_matx(n,2),position_matx(n,3),sigma_x, sigma_y);   % Gaussian wavepacket in the sublattice
                sign = symmetry*sign;                                                                       % Symmetric/ Antisymmetric wavepacket       
            % B sublattice
            else if H(n,n)<0 %&& position_matx(n,3)>y0-sqrt(3)
                position = [position_matx(n,2); position_matx(n,3)];                                        % Position of each state
                kick = exp(-1i*mom_kick*position)*exp(1i*mom_kick*[x0; y0]);                                % Momentum kick for the initial wavepcket
                Psi_0(n)=sign*kick*gauss_weight( position_matx(n,2),position_matx(n,3),sigma_x, sigma_y);   % Gaussian wavepacket in the sublattice
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
norm = sqrt(transpose(conj(Psi_0))*Psi_0);                                   % Normalisation
Psi_0 = Psi_0 / norm;                                                        % Normalisation
prob_density_0=abs(conj(Psi_0).*Psi_0);                                      % Initial probability density of the packet
for n=1:Ns_real                                                              % Initial overlap
    overlap_0(n)=abs(conj(transpose(eigenstates(:,n)))*Psi_0)^2;
end


%% Wave-packet evolution 
    
switch evolution
    
    case 'Exact'   
        time=0:T/3:t_exact;
        for j=1:length(time)
            % Initial step
            if j==1
                Psi_t(:,j)=Psi_0;
                prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));% Probability density of the packet
            end
            % Stepwise drive
            if step==1
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
            
            % Change of step
            step=step+1;
            if step>3
                step=1;
            end
        end
        
    case 'Stroboscopic'
        time=0:T:t_strob;
        for j=1:length(time)
            
            Psi_t(:,j)= Wavepacket_evolution( Psi_0, time(j), quasienergy, eigenstates);   % Evolution at time t, j many column vectors
            prob_density(:,j)=abs(conj(Psi_t(:,j)).*Psi_t(:,j));                           % Probability density of the packet
            
        end
        
    otherwise
        disp('Unknown evolution!')
end


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


% QUASIENERGY BANDS
figure(2) 
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
% print -depsc stripbands_[1,1].eps



% INITIAL OVERLAP
figure(3)
% Main figure
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
plot(x0, y0, '.k', 'MarkerSize', 16)  % Center of the wavepacket
xlim([-1 limx + 0.5])
switch edge
    case 'Upper'        
       ylim([limy-3 limy+0.5 ]) % upper edge
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
cont=1;
for n=1:Ns_real % We take the probability density of sites at the edge,
    switch edge
        case 'Upper'
            if position_matx(n,3)>y_sep*(Ny-2)-0.1
                evolution_matx(:,cont)=prob_density(n,:);
                cont=cont+1;
            end
        case 'Lower'
            if position_matx(n,3)<y_sep + 0.1
                evolution_matx(:,cont)=prob_density(n,:);
                cont=cont+1;
            end
            
        otherwise
            disp('Unknown edge!')
    end
end
  


figure(4)
s=surf(evolution_matx);
set(gcf,'Position',[900 400 600 500])
box on
view(2)
c=jet;
colormap(c(10:50,:));
caxis([0 0.1])
colorbar
xlabel('Edge Site','FontSize',25,'FontWeight','bold','Interpret','latex')
ylabel('Time [T]','FontSize',25,'FontWeight','bold','Interpret','latex')
h = colorbar;
set(get(h,'label'),'string','$\vert \Psi \vert^2$ [arbitrary units]','FontSize',20,'FontWeight','bold','Interpret','latex');
xlim([1 32])
ylim([1 length(time)])
s.EdgeColor = 'flat'
ax = gca;
ax.FontSize =25; 
ax.XTick = [1,10,20,30];




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
    set(gcf,'Position',[500 50 1100 700])  
    
    frame = getframe(gcf); % Gets the the plot as the frame for the video 
    writeVideo(vSim,frame);
    hold off;
    
    
    clf % Clear after each iteration
    
end 
close(vSim); % Close the animation


end






