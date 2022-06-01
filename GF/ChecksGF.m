%% CHECKS FOR THE GREEN'S FUNCTIONS


home
close all
clear all
addpath('functions')

%% Parameters

% Drive parameters
omega=4/3; %(Units of J)
T=2*pi/omega;
delta=0; %(Units of J)
lambda=3.4;
J_1=1;
J_2=1;
J_3=1;

% Time parameters
n_t=1000; % Resolution of the time vector used in the time evolution
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,3*T/3,n_t)]; % Time vector of the drive
%t=[linspace(0,T/4,n_t),linspace(T/4,2*T/4,n_t),linspace(2*T/4,3*T/4,n_t),linspace(4*T/4,T,n_t)]; % Time vector of the drive
dt=t(2)-t(1); % Differential time step

% Reciprocal Lattice Parameters (Check Sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% Full BZ (For Berry curvature)
kx=[-l:0.1:l];
ky=[-l:0.1:l];
dk=kx(2)-kx(1);

% Functions
deriv=@(A,B,h) [(A(1)-B(1))/(2*h) ; (A(2)-B(2))/(2*h)]; % Numerical vector derivative

% GF parameters
kx_GF=linspace(-pi, pi, 1000); ky_GF=0; dk_GF=kx_GF(2)-kx_GF(1); % Perpendicular momentum along x
eps_GF=linspace(-2*pi, 2*pi, 500); % Range of quasienergy BRANCHCUT 0 -PI--PI, BRANCHCUT 1 0--2PI
%HF=zeros(2,2,length(kx)); % Floquet hamiltonian definition
h0=zeros(1,length(kx)); % Floquet vector definition
hx=zeros(1,length(kx)); % Floquet vector definition
hy=zeros(1,length(kx)); % Floquet vector definition
hz=zeros(1,length(kx)); % Floquet vector definition

sigmax=[0, 1; 1, 0]; % Pauli matrix x
sigmay=[0, -s; s, 0]; % Pauli matrix y
sigmaz=[1, 0; 0, -1]; % Pauli matrix z

%% Berry curvature of the lower band
% To choose a smooth gauge we store all complex phases in the second entry
% of the Bloch eigenvector by multiplying with the conjugate of the first
% entry and then normalising.

% Calculation of the lower band eigenvectors
for j=1:length(kx);
    for z=1:length(ky)
        
        clear values vectors
        
        U=Time_evolution_bulk(delta, kx(j),ky(z), n_t, dt, J_1, J_2, J_3);
        HF=1i*logm(U); % Floquet Hamiltonian of the drive
        
        [A, B]=eig(HF); [values, vectors]=order_eigenvalues(B, A); % Ordering spectrum
        eps_1(j,z)=values(1); %  Select the lower band eigenvalue
        
        u_1(:,j,z)=vectors(:,1).*conj(vectors(1,1)); % By storing all complex phases in the second entry we fix the gauge
        u_1(:,j,z)=u_1(:,j,z)./sqrt(ctranspose(u_1(:,j,z))*u_1(:,j,z)); % Normalization
                       
    end
end

% Calculation of the Berry curvature
for j=1:length(kx)-1;
    for z=1:length(ky)-1;
        
        clear prod 
        
        % Derivatives of the Bloch state in the BZ
        der_x= deriv( u_1(:, j+1, z), u_1(:, j, z), dk);
        der_y= deriv( u_1(:, j, z+1), u_1(:, j, z), dk);
        prod= ctranspose(der_x)*der_y;
              
        Berry_curvature(j,z)=-2*imag( prod ); % Berry curvature at every point in the BZ
        
        % Avoiding singularities (Appearing bc of lack of grid definition)
        if Berry_curvature(j,z)>1.5 || Berry_curvature(j,z)<-0.8
            Berry_curvature(j,z)=NaN;
        end
        
        
    end
end


figure(1)
set(gcf,'Position',[250 300 1700 600])
subplot(2,4,[1,2,5,6])
s=surf(kx(1:end-1), ky(1:end-1), Berry_curvature);
box on
view(2)
h = colorbar;
set(get(h,'label'),'string','$\Omega(k_x, k_y)$','FontSize',15,'FontWeight','bold','Interpret','latex');
s.EdgeColor = 'flat'
xlabel('$k_x$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$k_y$','FontSize',15,'FontWeight','bold','Interpret','latex')
zlabel('$\Omega(k_x, k_y)$','FontSize',15,'FontWeight','bold','Interpret','latex')
xlim([-l l])
ylim([-l l])
title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda)])

subplot(2,4,[3,4,7,8])
s=surf(kx, ky, eps_1);
%set(gcf,'Position',[900 400 600 500])
box on
view(2)
s.EdgeColor = 'flat'
h = colorbar;
set(get(h,'label'),'string','$ET/J$','FontSize',15,'FontWeight','bold','Interpret','latex');
xlabel('$k_x$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$k_y$','FontSize',15,'FontWeight','bold','Interpret','latex')
zlabel('$ET/J$','FontSize',15,'FontWeight','bold','Interpret','latex')
xlim([-l l])
ylim([-l l])
title(['\omega/J=', num2str(omega),'    \Delta/J=', num2str(delta), '    \lambda=', num2str(lambda)])

%% Dummy hamiltonian check (MATRIX)
% We check that a topologically trivial hamiltonian with flatbands does
% indeed give no roots and the correct singularities in the GF by using
% direct matrix inverse.

for s=1:length(eps_GF)   
    
    % G(eps, kx, ky=0)
     for j=1:length(kx_GF)    
         
        clear M         
        HF=diag([pi/2, 3*pi/2]); % Trivial hamiltonian      
        M=(eps_GF(s)*eye(2)-HF); % Denom in GF
        F(:,:,j)=inv(M)./(2*pi); % GF for eps_GF(s) at kx_GF (/2pi for later integration)
        
     end
     
     % In gap GF G(eps, ky=0)
     sum=0; % Initial value of the integration
     for t=1:length(kx_GF)-1
     % Carefull as is a matrix integration we cannot just use sum!!!!!!!!!!
     sum=sum+(F(:,:,t)+F(:,:,t+1))./2*dk_GF; % Numerical integration of the matrix F for -pi<kx<pi
        
     end

    % GF eigenvalues
    eigenvalues(:,s)=eig(sum); % Eigenvalues of the Green's function
       
end


figure(3)
hold on
box on
plot(eps_GF, eigenvalues(1,:),'.r')
plot(eps_GF, eigenvalues(2,:),'.b')
plot(ones(200).*pi/2,linspace(-10,10,200)) % Correct divergences
plot(ones(200).*3*pi/2,linspace(-10,10,200)) % Correct divergences
plot(eps_GF, 0*ones(1, length(eps_GF)),'r') %, '    $\lambda$=', num2str(lambda)
title('Greens function eigenvalues ','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
xlim([eps_GF(1) eps_GF(end)])
ax = gca;
ax.FontSize =15; 

%% Dummy hamiltonian check (PAULI)
% We check that a topologically trivial hamiltonian with flatbands does
% indeed give no roots and the correct singularities in the GF by using
% the Puali decomposition a=eps-h0 b=h*sigma.

% Hamiltonian decomposition
for i=1:length(kx_GF)
    
    HF=diag([pi/2, 3*pi/2]); % Trivial hamiltonian
    
    % Pauli decomposition
    hx(i)=real(HF(2,1)); h0(i)=(HF(1,1)+HF(2,2))/2;
    hy(i)=imag(HF(2,1)); hz(i)=(HF(1,1)-HF(2,2))/2;
    h2(i)=hx(i)^2+hy(i)^2+hz(i)^2;
    
end

%Calculation of the Green's function
for s=1:length(eps_GF)
    
    % Integrands for g, gx, gy, gz
    for j=1:length(kx_GF)
    fx(j)=hx(j)/(((eps_GF(s)-h0(j))^2-h2(j))*2*pi);
    fy(j)=hy(j)/(((eps_GF(s)-h0(j))^2-h2(j))*2*pi);
    fz(j)=hz(j)/(((eps_GF(s)-h0(j))^2-h2(j))*2*pi);
    f(j)=(eps_GF(s)-h0(j))/(((eps_GF(s)-h0(j))^2-h2(j))*2*pi);
    end
    
    % Integrals for g, g0, gx, gy, gz
    Gx(s)=sum(fx(1:length(kx_GF)-1)+ fx(2:length(kx_GF)))/2*dk_GF;
    Gy(s)=sum(fy(1:length(kx_GF)-1)+ fy(2:length(kx_GF)))/2*dk_GF;
    Gz(s)=sum(fz(1:length(kx_GF)-1)+ fz(2:length(kx_GF)))/2*dk_GF;
    G(s)=sum(f(1:length(kx_GF)-1)+ f(2:length(kx_GF)))/2*dk_GF;
    
    % Complete Green's function
    G_tot=G(s)*eye(2)+Gx(s)*sigmax+Gy(s)*sigmay+Gz(s)*sigmaz;
    eigenvalues(:,s)=eig(G_tot); % Eigenvalues of the Green's function
    
end


figure(4)
hold on
box on
plot(eps_GF, eigenvalues(1,:),'.b')
plot(eps_GF, eigenvalues(2,:),'.b')
plot(ones(200).*pi/2,linspace(-10,10,200)) % Correct divergences
plot(ones(200).*3*pi/2,linspace(-10,10,200))% Correct divergences
plot(eps_GF, 0*ones(1, length(eps_GF)),'r') %, '    $\lambda$=', num2str(lambda)
title('Greens function eigenvalues ','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
xlim([eps_GF(1) eps_GF(end)])
ax = gca;
ax.FontSize =15; 

%% Check Double


for s=1:length(eps_GF)   
    s
    % G(eps, kx, ky=0)
     for j=1:length(kx_GF)    
         
         clear HF values vectors A B M
         
      H=diag([pi/2, 3*pi/2]);
      
      HF=[zeros(2), expm(-i*H);
          ctranspose(expm(-i*H)), zeros(2)];
         
         M=(eps_GF(s)*eye(4)-HF); % Denom in GF
         F(:,:,j)=inv(M)./(2*pi); % GF for eps_GF(s) at kx_GF (/2pi for later integration)
         
     end
     
     % In gap GF G(eps, ky=0)
     sum=0; % Initial value of the integration
     for t=1:length(kx_GF)-1
     % Careful as is a matrix integration we cannot just use sum!!!!!!!!!!
     sum=sum+(F(:,:,t)+F(:,:,t+1))./2*dk_GF; % Numerical integration of the matrix F for -pi<kx<pi
        
     end

    % GF eigenvalues
    eigenvalues(:,s)=eig(sum); % Eigenvalues of the Green's function
       
end



figure(2)
hold on
box on
plot(eps_GF, eigenvalues(1,:),'.b','MarkerSize',15)
plot(eps_GF, eigenvalues(2,:),'.b','MarkerSize',15)
plot(eps_GF, eigenvalues(3,:),'.b','MarkerSize',15)
plot(eps_GF, eigenvalues(4,:),'.b','MarkerSize',15)
plot(eps_GF, 1./(eps_GF-1), 'c', eps_GF, 1./(eps_GF+1), 'k')
plot(eps_GF, 0*ones(1, length(eps_GF)),'r') %, '    $\lambda$=', num2str(lambda)
title('Greens function eigenvalues ','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$ET$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
xlim([eps_GF(1) eps_GF(end)])
ax = gca;
ax.FontSize =15; 




