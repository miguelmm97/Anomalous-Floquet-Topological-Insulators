%% 3-STEP PERIODIC DRIVING
clear all
home

% Trial Basis
u1_0=[0,1];
u2_0=[-1,0];
basis_0=[u1_0', u2_0'];
Periods=[2*pi/3];

opts = odeset('RelTol',1e-2,'AbsTol',1e-2);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
%% Gamma --> K

kx=[0:0.1:l];
ky=zeros(length(kx));

for j=1:length(Periods)
T=Periods(j);
% Driving step 1
for i=1:length(kx);
    
    H=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    
    [t,y1]=ode45(f1,[0,T/3],u1_0,opts);
    [t,y2]=ode45(f2,[0,T/3],u2_0,opts);
    
     % Resulting basis after 1st driving step
     u1_1(i,:)=y1(end,:);
     u2_1(i,:)=y2(end,:);
end

% Driving step 2
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[T/3,2*T/3],u1_1(i,:),opts);
    [t,y2]=ode45(f2,[T/3,2*T/3],u2_1(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_2(i,:)=y1(end,:);
     u2_2(i,:)=y2(end,:);
end

% Driving step 3
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[2*T/3,T],u1_2(i,:),opts);
    [t,y2]=ode45(f2,[2*T/3,T],u2_2(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_3(i,:)=y1(end,:);
     u2_3(i,:)=y2(end,:);
          
     % Time evolution operator
   
     U=zeros(2,2);
     U(1,1)=conj(u1_0)*u1_3(i,:)';
     U(2,1)=conj(u2_0)*u1_3(i,:)';
     U(1,2)=conj(u1_0)*u2_3(i,:)';
     U(2,2)=conj(u2_0)*u2_3(i,:)';
     
     
     %U=conj(basis_0)'*[u1_3(i,:)', u2_3(i,:)'];
     U_diag=eig(U);
     % Quasi-Energy Bands
     eps_Gamma_K_1(i,j)=1i*log(U_diag(1))./Periods(j);
     eps_Gamma_K_2(i,j)=1i*log(U_diag(2))./Periods(j);
     
end

end

%% K -->M


kx=[l:-0.01:a];
ky=tan(120*pi/180)*kx-(l*tan(120*pi/180));

for j=1:length(Periods)
T=Periods(j);
% Driving step 1
for i=1:length(kx);
    
    H=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[0,T/3],u1_0,opts);
    [t,y2]=ode45(f2,[0,T/3],u2_0,opts);
    
     % Resulting basis after 1st driving step
     u1_1(i,:)=y1(end,:);
     u2_1(i,:)=y2(end,:);
end

% Driving step 2
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[T/3,2*T/3],u1_1(i,:),opts);
    [t,y2]=ode45(f2,[T/3,2*T/3],u2_1(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_2(i,:)=y1(end,:);
     u2_2(i,:)=y2(end,:);
end

% Driving step 3
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[2*T/3,T],u1_2(i,:),opts);
    [t,y2]=ode45(f2,[2*T/3,T],u2_2(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_3(i,:)=y1(end,:);
     u2_3(i,:)=y2(end,:);
          
     % Time evolution operator
     U=zeros(2,2);
     U(1,1)=conj(u1_0)*u1_3(i,:)';
     U(2,1)=conj(u2_0)*u1_3(i,:)';
     U(1,2)=conj(u1_0)*u2_3(i,:)';
     U(2,2)=conj(u2_0)*u2_3(i,:)';
     
     U_diag=eig(U);
     % Quasi-Energy Bands
     eps_K_M_1(i,j)=1i*log(U_diag(1))./Periods(j);
     eps_K_M_2(i,j)=1i*log(U_diag(2))./Periods(j);
end
end
%% M --> Gamma


kx=[a:-0.01:0];
ky=tan(30*pi/180)*kx;

for j=1:length(Periods)
    T=Periods(j);
% Driving step 1
for i=1:length(kx);
    
    H=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[0,T/3],u1_0,opts);
    [t,y2]=ode45(f2,[0,T/3],u2_0,opts);
    
     % Resulting basis after 1st driving step
     u1_1(i,:)=y1(end,:);
     u2_1(i,:)=y2(end,:);
end

% Driving step 2
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[T/3,2*T/3],u1_1(i,:),opts);
    [t,y2]=ode45(f2,[T/3,2*T/3],u2_1(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_2(i,:)=y1(end,:);
     u2_2(i,:)=y2(end,:);
end

% Driving step 3
for i=1:length(kx);
    
    H=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    
    % Schrodinger's Equation for the 1st driving stage
    f1= @(t,y1) [-1i*(H(1,1)*y1(1)+H(1,2)*y1(2)); -1i*(H(2,1)*y1(1)+H(2,2)*y1(2))];
    f2= @(t,y2) [-1i*(H(1,1)*y2(1)+H(1,2)*y2(2)); -1i*(H(2,1)*y2(1)+H(2,2)*y2(2))];
    
    [t,y1]=ode45(f1,[2*T/3,T],u1_2(i,:),opts);
    [t,y2]=ode45(f2,[2*T/3,T],u2_2(i,:),opts);
    
     % Resulting basis after 1st driving step
     u1_3(i,:)=y1(end,:);
     u2_3(i,:)=y2(end,:);
          
     % Time evolution operator
     U=zeros(2,2);
     U(1,1)=conj(u1_0)*u1_3(i,:)';
     U(2,1)=conj(u2_0)*u1_3(i,:)';
     U(1,2)=conj(u1_0)*u2_3(i,:)';
     U(2,2)=conj(u2_0)*u2_3(i,:)';
     
     U_diag=eig(U);
     % Quasi-Energy Bands
     eps_M_Gamma_1(i,j)=1i*log(U_diag(1))./Periods(j);
     eps_M_Gamma_2(i,j)=1i*log(U_diag(2))./Periods(j);
end
end
%% QUASI-ENERGY BAND STRUCTURE

K=l;
M=sqrt(b^2+(l-a)^2);
Gamma=2*pi/3;

k1=linspace(0,K,length(eps_Gamma_K_1));
k2=linspace(K,K+M,length(eps_K_M_1));
k3=linspace(K+M,K+M+Gamma,length(eps_M_Gamma_1));


%colours=['.b','.r','.g','.c'];
close all
figure(1)
box on
hold on
for j=1:length(Periods)
    if j==1
        c='.b';
    else if j==2
            c='.r';
        else if j==3
                c='.g';
            else
                c='.k';
            end
        end
    end
plot(k1,eps_Gamma_K_1(:,j).*Periods(j),c,k1,eps_Gamma_K_2(:,j).*Periods(j),c)
plot(k2,eps_K_M_1(:,j).*Periods(j),c,k2,eps_K_M_2(:,j).*Periods(j),c)
plot(k3,eps_M_Gamma_1(:,j).*Periods(j),c',k3,eps_M_Gamma_2(:,j).*Periods(j),c)
end
plot(K*ones(length([-3:3])),[-3:3],'--k')
plot(K+M*ones(length([-3:3])),[-3:3],'--k')
hold off

%legend('Numerical dispersion')

ax=gca;
ax.XTick=[0,K,K+M,K+M+Gamma];
ax.XTickLabel={'\Gamma','K','M','\Gamma'}
ax = gca;
ax.FontSize = 16;
xlim([0 K+M+Gamma])
ylim([-pi pi])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$ET/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off




