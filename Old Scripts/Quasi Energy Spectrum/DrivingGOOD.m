%% 3-STEP PERIODIC DRIVING
clear all
home



Periods=[3*pi/3];
n_t=1000;
T=Periods(1);
t=[linspace(0,T/3,n_t),linspace(T/3,2*T/3,n_t),linspace(2*T/3,T,n_t)];
dt=t(2)-t(1);

% Reciprocal Lattice Parameters
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
%% Gamma --> K

kx=[0:0.01:l];
ky=zeros(length(kx),1);

for i=1:length(kx);
    
    H1=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];

    
    H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
   
    
    H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
      
     U1=(expm(-1i*dt*H1))^n_t;
     U2=(expm(-1i*dt*H2))^n_t;
     U3=(expm(-1i*dt*H3))^n_t;
     prod=U1*U2*U3;
     
     U_diag=eig(prod);
     % Quasi-Energy Bands
     eps_Gamma_K_1(i)=1i*log(U_diag(1))./T;
     eps_Gamma_K_2(i)=1i*log(U_diag(2))./T;
     
end

%% K -->M

kx=[l:-0.01:a];
ky=tan(120*pi/180)*kx-(l*tan(120*pi/180));

for i=1:length(kx);
    
    H1=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];

    
    H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
   
    
    H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
      
     U1=(expm(-1i*dt*H1))^n_t;
     U2=(expm(-1i*dt*H2))^n_t;
     U3=(expm(-1i*dt*H3))^n_t;
     prod=U1*U2*U3;
     
     U_diag=eig(prod);
     % Quasi-Energy Bands
     eps_K_M_1(i)=1i*log(U_diag(1))./T;
     eps_K_M_2(i)=1i*log(U_diag(2))./T;
     
end

    

%% M --> Gamma

kx=[a:-0.01:0];
ky=tan(30*pi/180)*kx;

for i=1:length(kx);
    
    H1=[0, exp(1i*ky(i));
        exp(-1i*ky(i)), 0];

    
    H2=[0, exp(-1i*0.5*(sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(sqrt(3)*kx(i)+ky(i))), 0];
   
    
    H3=[0, exp(-1i*0.5*(-sqrt(3)*kx(i)+ky(i)));
        exp(1i*0.5*(-sqrt(3)*kx(i)+ky(i))), 0];
      
     U1=(expm(-1i*dt*H1))^n_t;
     U2=(expm(-1i*dt*H2))^n_t;
     U3=(expm(-1i*dt*H3))^n_t;
     prod=U1*U2*U3;
     
     U_diag=eig(prod);
     % Quasi-Energy Bands
     eps_M_Gamma_1(i)=1i*log(U_diag(1))./T;
     eps_M_Gamma_2(i)=1i*log(U_diag(2))./T;
     
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
plot(k1,eps_Gamma_K_1*T,c,k1,eps_Gamma_K_2*T,c)
plot(k2,eps_K_M_1*T,c,k2,eps_K_M_2*T,c)
plot(k3,eps_M_Gamma_1*T,c',k3,eps_M_Gamma_2*T,c)
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

