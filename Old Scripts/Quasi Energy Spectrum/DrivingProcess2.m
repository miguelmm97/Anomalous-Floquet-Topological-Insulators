%% 3-STEP PERIODIC DRIVING
clear all
home

% Trial Basis
u1_0=[1;0];
u2_0=[0;1];
basis_0=[u1_0, u2_0];

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;

%% Gamma-->K

kx=[0:0.01:1];
ky=zeros(length(kx));

% 1st Driving Step
for i=1:length(kx)
    
H1=[0, exp(1i*ky(i));
   exp(-1i*ky(i)), 0];

[V1,D1]=eig(H1);

E1=diag(D1);

u1_1(:,i)=[(conj(V1(:,1))'*u1_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u1_0)*exp(-1i*E1(2)/3)];
u2_1(:,i)=[(conj(V1(:,1))'*u2_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u2_0)*exp(-1i*E1(2)/3)];
end


% 2nd Driving Step

for i=1:length(kx)

H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
  exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    

[V2,D2]=eig(H2);

E2=diag(D2);

phi1_2=inv(V2)*V1(:,1);
phi2_2=inv(V2)*V1(:,2);

u1_1_2=u1_1(1,i)*phi1_2+u1_1(2,i)*phi2_2;
u2_1_2=u2_1(1,i)*phi1_2+u2_1(2,i)*phi2_2;

u1_2(:,i)=[(conj(V2(:,1))'*u1_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u1_1_2)*exp(-1i*2*E2(2)/3)];
u2_2(:,i)=[(conj(V2(:,1))'*u2_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u2_1_2)*exp(-1i*2*E2(2)/3)];
  
end

% 3rd Driving step

for i=1:length(kx)
    
H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
    exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    

[V3,D3]=eig(H3);

E3=diag(D3);

phi1_3=inv(V3)*V2(:,1);
phi2_3=inv(V3)*V2(:,2);

u1_2_3=u1_2(1,i)*phi1_3+u1_2(2,i)*phi2_3;
u2_2_3=u2_2(1,i)*phi1_3+u2_2(2,i)*phi2_3;

u1_3(:,i)=[(conj(V3(:,1))'*u1_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u1_2_3)*exp(-1i*E3(2))];
u2_3(:,i)=[(conj(V3(:,1))'*u2_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u2_2_3)*exp(-1i*E3(2))];

  
end

% Time Evolution operator
u1_0_3=inv(V3)*u1_0;
u2_0_3=inv(V3)*u2_0;

for i=1:length(kx)
    
U=zeros(2,2);
U(1,1)=conj(u1_0_3)'*u1_3(:,i);
U(2,1)=conj(u2_0_3)'*u1_3(:,i);
U(1,2)=conj(u1_0_3)'*u2_3(:,i);
U(2,2)=conj(u2_0_3)'*u2_3(:,i);

U_diag=eig(U);
     
% Quasi-Energy Bands
eps_Gamma_K_1(i)=1i*log(U_diag(1));
eps_Gamma_K_2(i)=1i*log(U_diag(2));    
end


%% K -->M

kx=[l:-0.01:a];
ky=tan(120*pi/180)*kx-(l*tan(120*pi/180));

% 1st Driving Step
for i=1:length(kx)
    
H1=[0, exp(1i*ky(i));
   exp(-1i*ky(i)), 0];

[V1,D1]=eig(H1);

E1=diag(D1);

u1_1(:,i)=[(conj(V1(:,1))'*u1_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u1_0)*exp(-1i*E1(2)/3)];
u2_1(:,i)=[(conj(V1(:,1))'*u2_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u2_0)*exp(-1i*E1(2)/3)];
end

% 2nd Driving Step

for i=1:length(kx)

H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
  exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    

[V2,D2]=eig(H2);

E2=diag(D2);

phi1_2=inv(V2)*V1(:,1);
phi2_2=inv(V2)*V1(:,2);

u1_1_2=u1_1(1,i)*phi1_2+u1_1(2,i)*phi2_2;
u2_1_2=u2_1(1,i)*phi1_2+u2_1(2,i)*phi2_2;

u1_2(:,i)=[(conj(V2(:,1))'*u1_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u1_1_2)*exp(-1i*2*E2(2)/3)];
u2_2(:,i)=[(conj(V2(:,1))'*u2_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u2_1_2)*exp(-1i*2*E2(2)/3)];
  
end

% 3rd Driving step

for i=1:length(kx)
    
H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
    exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    

[V3,D3]=eig(H3);

E3=diag(D3)

phi1_3=inv(V3)*V2(:,1);
phi2_3=inv(V3)*V2(:,2);

u1_2_3=u1_2(1,i)*phi1_3+u1_2(2,i)*phi2_3;
u2_2_3=u2_2(1,i)*phi1_3+u2_2(2,i)*phi2_3;

u1_3(:,i)=[(conj(V3(:,1))'*u1_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u1_2_3)*exp(-1i*E3(2))];
u2_3(:,i)=[(conj(V3(:,1))'*u2_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u2_2_3)*exp(-1i*E3(2))];

  
end

% Time Evolution operator
u1_0_3=inv(V3)*u1_0;
u2_0_3=inv(V3)*u2_0;

for i=1:length(kx)
    
U=zeros(2,2);
U(1,1)=conj(u1_0_3)'*u1_3(:,i);
U(2,1)=conj(u2_0_3)'*u1_3(:,i);
U(1,2)=conj(u1_0_3)'*u2_3(:,i);
U(2,2)=conj(u2_0_3)'*u2_3(:,i);

U_diag=eig(U);
     
% Quasi-Energy Bands
eps_K_M_1(i)=1i*log(U_diag(1));
eps_K_M_2(i)=1i*log(U_diag(2));    
end


%% M --> Gamma

kx=[a:-0.01:0];
ky=tan(30*pi/180)*kx;



% 1st Driving Step
for i=1:length(kx)
    
H1=[0, exp(1i*ky(i));
   exp(-1i*ky(i)), 0];

[V1,D1]=eig(H1);

E1=diag(D1);

u1_1(:,i)=[(conj(V1(:,1))'*u1_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u1_0)*exp(-1i*E1(2)/3)];
u2_1(:,i)=[(conj(V1(:,1))'*u2_0)*exp(-1i*E1(1)/3);(conj(V1(:,2))'*u2_0)*exp(-1i*E1(2)/3)];
end

% 2nd Driving Step

for i=1:length(kx)

H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
  exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
    

[V2,D2]=eig(H2);

E2=diag(D2);

phi1_2=inv(V2)*V1(:,1);
phi2_2=inv(V2)*V1(:,2);

u1_1_2=u1_1(1,i)*phi1_2+u1_1(2,i)*phi2_2;
u2_1_2=u2_1(1,i)*phi1_2+u2_1(2,i)*phi2_2;

u1_2(:,i)=[(conj(V2(:,1))'*u1_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u1_1_2)*exp(-1i*2*E2(2)/3)];
u2_2(:,i)=[(conj(V2(:,1))'*u2_1_2)*exp(-1i*2*E2(1)/3);(conj(V2(:,2))'*u2_1_2)*exp(-1i*2*E2(2)/3)];
  
end

% 3rd Driving step

for i=1:length(kx)
    
H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
    exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
    

[V3,D3]=eig(H3);

E3=diag(D3);

phi1_3=inv(V3)*V2(:,1);
phi2_3=inv(V3)*V2(:,2);

u1_2_3=u1_2(1,i)*phi1_3+u1_2(2,i)*phi2_3;
u2_2_3=u2_2(1,i)*phi1_3+u2_2(2,i)*phi2_3;

u1_3(:,i)=[(conj(V3(:,1))'*u1_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u1_2_3)*exp(-1i*E3(2))];
u2_3(:,i)=[(conj(V3(:,1))'*u2_2_3)*exp(-1i*E3(1));(conj(V3(:,2))'*u2_2_3)*exp(-1i*E3(2))];

  
end

% Time Evolution operator
u1_0_3=inv(V3)*u1_0;
u2_0_3=inv(V3)*u2_0;

for i=1:length(kx)
    
U=zeros(2,2);
U(1,1)=conj(u1_0_3)'*u1_3(:,i);
U(2,1)=conj(u2_0_3)'*u1_3(:,i);
U(1,2)=conj(u1_0_3)'*u2_3(:,i);
U(2,2)=conj(u2_0_3)'*u2_3(:,i);

U_diag=eig(U);
     
% Quasi-Energy Bands
eps_M_Gamma_1(i)=1i*log(U_diag(1));
eps_M_Gamma_2(i)=1i*log(U_diag(2));    
end


%% QUASI-ENERGY BAND STRUCTURE

K=l;
M=sqrt(b^2+(l-a)^2);
Gamma=2*pi/3;

k1=linspace(0,K,length(eps_Gamma_K_1));
k2=linspace(K,K+M,length(eps_K_M_1));
k3=linspace(K+M,K+M+Gamma,length(eps_M_Gamma_1));

close all
figure(1)
box on
hold on
plot(k1,eps_Gamma_K_1,'-r',k1,eps_Gamma_K_2,'-b')
plot(k2,eps_K_M_1,'-r',k2,eps_K_M_2,'-b')
plot(k3,eps_M_Gamma_1,'-r',k3,eps_M_Gamma_2,'-b')
%plot(zeros(length([0:3])),[0:3],'--k')

%legend('Numerical dispersion')

ax=gca;
ax.XTick=[0,K,K+M,K+M+Gamma];
ax.XTickLabel={'\Gamma','K','M','\Gamma'}
ax = gca;
ax.FontSize = 16;
xlim([0 K+M+Gamma])
xlabel('$ka$','FontSize',18,'FontWeight','bold','Interpret','latex')
ylabel('$E/J$','FontSize',18,'FontWeight','bold','Interpret','latex')
hold off

