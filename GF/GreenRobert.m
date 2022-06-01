%% GREEN'S FUNCTIONS ROBERT PAPER (CHECK OF WORKING)

% In this code we compute the in-gap Green's functions for the Floquet topological
% phases following the thoery described by Robert.


home
close all
clear all
addpath('functions')

%% Parameters

M=1;
B=1;

% Reciprocal Lattice Parameters (Check Sketch in the project notes)
l=sqrt((4*(pi^2)/9)/(3/4));
a=2*pi*cos(30*pi/180)/3;
b=2*pi*sin(30*pi/180)/3;
s=l*cos(60*pi/180);
f=l*sin(60*pi/180);

% Definitions
kx=linspace(-pi, pi, 2000); ky=0; dk=kx(2)-kx(1);
eps=linspace(-pi, pi, 2000); 
HF=zeros(2,2,length(kx));
h0=zeros(1,length(kx));
hx=zeros(1,length(kx));
hy=zeros(1,length(kx));
hz=zeros(1,length(kx));

sigma0=[1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1];
sigmax=[0, 1, 0, 0; 1, 0, 0, 0; 0, 0, 0, 1; 1, 0, 0, 0];
sigmay=[0, -i, 0, 0; i, 0, 0, 0; 0, 0, 0, 1; 1, 0, 0, 0];


h0=M-2*B*(2-cos(kx)-cos(ky));
hx=sin(kx);
hy=sin(ky);
h2=+hx.^2+hy.^2+h0.^2

 

%% Calculation of the Green's function



for i=1:length(eps)
    f0=h0./((eps(i)^2-h2)*2*pi);
    fx=hx./((eps(i)^2-h2)*2*pi);
    fy=hy./((eps(i)^2-h2)*2*pi);
    f=eps(i)./((eps(i)^2-h2)*2*pi); 
       
    G0(i)=sum(f0(1:length(kx)-1)+ f0(2:length(kx)))/2*dk;
    Gx(i)=sum(fx(1:length(kx)-1)+ fx(2:length(kx)))/2*dk;
    Gy(i)=sum(fy(1:length(kx)-1)+ fy(2:length(kx)))/2*dk;
    G(i)=sum(f(1:length(kx)-1)+ f(2:length(kx)))/2*dk;
     
    G_tot(:,:,i)=G(i)*eye(4)+G0(i)*sigma0+Gx(i)*sigmax+Gy(i)*sigmay;    
    eigenvalues(:,i)=eig(G_tot(:,:,i));
end

figure(1)
hold on
box on
plot(eps, eigenvalues(1,:),'.b')
plot(eps, eigenvalues(4,:),'.b')
plot(eps, -G0+G,'or')
plot(eps, G0+G,'or')
plot(eps, 0*ones(1, length(eps)),'r')
%plot(0*ones(1, length(eps)), linspace(-500,500,length(eps)),'r')

title('In-gap Greens functions eigenvalues','FontSize',15,'FontWeight','bold','Interpret','latex')
xlabel('$\epsilon$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylabel('$\lambda$','FontSize',15,'FontWeight','bold','Interpret','latex')
ylim([-2 2])
%xlim([eps(1) eps(end)])
xlim([-0.8 0.8])
ax = gca;
ax.FontSize =15; 
%ax.XTick = [-pi/T, 0, pi/T];
%ax.XTickLabel = {'-\pi/T','0','\pi/T'};
% ax.YTick = [-pi, -2, -1, 0, 1, 2, pi];
% ax.YTickLabel = {'-\pi','-2','-1','0','1','2','\pi'};


