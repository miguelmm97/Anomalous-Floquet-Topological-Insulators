% Here I calculate the Hopf invariant of the state following a quench    
% applying it to honeycomb, 3stage step-wise driving  

% I calculate the Hopf invariant over N number of periods, allowing for doubling yani  
% I also apply Return Map after Nperiods

% I don't go to Bloch sphere, but instead consider the evolution in the sublattice basis  

% initial condition z_i=[1; 0],  I analytically write the rotation matrices
% U=cos(alpha/2) -i*sin(alpha/2)*(h),  where h=2x2 matrix, alpha=2*|h(k)|*t   ***************** be careful,  delta~=0 ise |h(k)|~=1,  so alpha=2t olmicak!!!!!!  
% z(kx,ky,t)= U*z_i,   so for every time stage I provide this z as function_handle 

function StepWiseHxgn_HopfOfState_ReturnMap
global J1 J2 J3 T deltaA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- INPUT --------------------------------
J=1;     % everything is in terms of J
w=1.4;  % driving frequency omega, hbar*omega/J, T=2pi/w 
delta=0;   % detuning. for initial sublattice bias, divided equally btwn A and B sublattices so delta_B=delta/2, delta_A=-delta_B

% N=1;    % your period will be N*T

J1=1;   % if you want anisotropy along the three hoppings
J2=1;
J3=3;

dk=.03;     % increment in dkx dky for integral
tN=20;    % how many steps I divide a single time stage, which is T/3

derk=.001;   % this is to take derivative in kx, ky
dert=.0005;   %for time derivatives
% branchcut=0;   %I'll shift the branchcut of the log, 0=>Hopf of zone edge, 1=>hopf of center
% -----------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kxRange=0:dk:2*pi;     % BZ is a rectangle!, so no Jacobian needed, BE CAREFUL!!
kyRange=0:dk:(4*pi/sqrt(3));  
kNx=length(kxRange);  kNy=length(kyRange);

epsilon=1e-5;
kxRange=linspace(0+epsilon,2*pi-epsilon,kNx);    dx=kxRange(2)-kxRange(1);   % here I include both 0 and 2pi in the Range up to epsilon, in case it causes any problems...
kyRange=linspace(0+epsilon,4*pi/sqrt(3)-epsilon,kNy);   dy=kyRange(2)-kyRange(1);



T=2*pi/w;   % period of the drive
deltaB=delta/2;  deltaA=-deltaB;  %putting offset on equal footing to sublatices

dt=T/3/tN;   % time steps for time integral
if dert>dt,  warning('derivative is bigger than the integral range!!! wanna continue stupidity?'); pause,  end




hmag1=sqrt(deltaB^2+J1^2);  % and h1=h2=h3 amplitudes for the rotations
hmag2=sqrt(deltaB^2+J2^2);
hmag3=sqrt(deltaB^2+J3^2);

U1=@(tt)       (cos(hmag1*tt)*eye(2)-1i*sin(hmag1*tt)* [deltaA -J1; -J1 deltaB]/hmag1);  
U2=@(kx,ky,tt) (cos(hmag2*tt)*eye(2)-1i*sin(hmag2*tt)* [deltaA -J2*exp(1i/2*(kx+sqrt(3)*ky)); -J2*exp(-1i/2*(kx+sqrt(3)*ky)) deltaB]/hmag2);   % BE CAREFUL when calling, THESE t ARE THE EVOLUTION TIMES, delta_t YANI, o stage icinde ne kadar sureyle evolve ettigin
U3=@(kx,ky,tt) (cos(hmag3*tt)*eye(2)-1i*sin(hmag3*tt)* [deltaA -J3*exp(1i/2*(-kx+sqrt(3)*ky)); -J3*exp(-1i/2*(-kx+sqrt(3)*ky)) deltaB]/hmag3);



z_i=[1; 0];    %initial state, in sublattice basis [1; 0]

z_fun=@(kx,ky,t) U1(t)*z_i;   % wave function during the first stage. I need to call this function;  %during the 1stage, I'll calculate the w.f.s by using this function

tic
%% ------------------- 1st period 0-->T ---------------------------------
% filename=[pwd '/data/Hopf_w' num2str(w) '_d' num2str(delta) '_dk' num2str(dk) '_tN' num2str(tN) '_derk' num2str(derk) '_dert' num2str(dert) '_J1' num2str(J1) '_J2' num2str(J2) '_J3' num2str(J3) '.mat'];
% if exist(filename,'file')==2    % if the iterated wf already exists, load
%     disp('Hopf1 exists, loading...')
%     load(filename);
% else
    disp('Hopf1 doesn''t exist, calculating...')

    intgd_t=zeros(1,3*tN+1);   % data points at each time step, including time=0 and T    
    for cnt=0:3*tN  
        switch cnt
            case tN  %switch to the 2nd stage
                z_fun=@(kx,ky,t) U2(kx,ky,t-T/3)* U1(T/3) *z_i;     %function during the 2nd stage. t=how much you evolve only during the 2nd stage
%                 cnt, %pause
            case 2*tN  %3rd stage
                z_fun=@(kx,ky,t) U3(kx,ky,t-2*T/3)* U2(kx,ky,T/3)* U1(T/3) *z_i;     
%                 cnt, %pause
        end   %now I have the wave functions defined in z_fun(kx,ky,t) for each stage




        intgd_y=zeros(1,kNy);  % data points in dky integral
        for cny=1:kNy   % to calculate integrand for dky integral   
            y=kyRange(cny);

            intgd_x=zeros(1,kNx);  % data points in dkx integral
            for cnx=1:kNx    % to calculate integrand of dkx integral   
                x=kxRange(cnx);             

                z0=z_fun(x,y,cnt*dt);
                dxz=( z_fun(x+derk/2,y,cnt*dt)-z_fun(x-derk/2,y,cnt*dt) )/derk;    
                dyz=( z_fun(x,y+derk/2,cnt*dt)-z_fun(x,y-derk/2,cnt*dt) )/derk;
                dtz=( z_fun(x,y,cnt*dt+dert/2)-z_fun(x,y,cnt*dt-dert/2) )/dert;                   
                if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end

                % --------------- alternative for the levi-civita sum
                intgd_x(cnx)= 2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));     
            end      

            intgd_y(cny)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral            
        end       

        intgd_t(cnt+1)=sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy;   % here I took dky integral,  index t+1 bcs loop count starts from zero  



        % this part might be commented out maybe-discontinuity may not matter much for most cases....              
        % -------- now I'm taking care of the discontinuity btwn stages. I'll calculate the dt integral at the end of each stage and save in-----------   
        if cnt==(tN-1)  %so I need to calculate the last data point of 1st stage before changing the Hamiltonian, z hala 1st stage'in w.f.'i
            intgd_y=zeros(1,kNy);  % data points in dky integral
            for cny=1:kNy   % to calculate integrand for dky integral   
                y=kyRange(cny);
                intgd_x=zeros(1,kNx);  % data points in dkx integral
                for cnx=1:kNx    % to calculate integrand of dkx integral   
                    x=kxRange(cnx);  

                    z0=z_fun(x,y,(cnt+1)*dt);     % I'm passing (t+1)*dt here
                    dxz=( z_fun(x+derk/2,y,(cnt+1)*dt)-z_fun(x-derk/2,y,(cnt+1)*dt) )/derk;
                    dyz=( z_fun(x,y+derk/2,(cnt+1)*dt)-z_fun(x,y-derk/2,(cnt+1)*dt) )/derk;
                    dtz=( z_fun(x,y,(cnt+1)*dt+dert/2)-z_fun(x,y,(cnt+1)*dt-dert/2) )/dert;                
                    if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end

                    intgd_x(cnx)= 2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));                
                end

                intgd_y(cny)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral    
            end
            discont_t(1)=sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy;   % here I took dky integral, for data point at t=T/4 by using H1

        elseif cnt==(2*tN-1)   %last data point of 2nd stg
            intgd_y=zeros(1,kNy);  % data points in dky integral
            for cny=1:kNy   % to calculate integrand for dky integral    
                y=kyRange(cny);
                intgd_x=zeros(1,kNx);  % data points in dkx integral
                for cnx=1:kNx    % to calculate integrand of dkx integral   
                    x=kxRange(cnx); 

                    z0=z_fun(x,y,(cnt+1)*dt);
                    dxz=( z_fun(x+derk/2,y,(cnt+1)*dt)-z_fun(x-derk/2,y,(cnt+1)*dt) )/derk;
                    dyz=( z_fun(x,y+derk/2,(cnt+1)*dt)-z_fun(x,y-derk/2,(cnt+1)*dt) )/derk;
                    dtz=( z_fun(x,y,(cnt+1)*dt+dert/2)-z_fun(x,y,(cnt+1)*dt-dert/2) )/dert;                
                    if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end

                    intgd_x(cnx)= 2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz)); 
                end

                intgd_y(cny)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral    
            end
            discont_t(2)=sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy;   % here I took dky integral, for data point at t=2*T/4 by using H2

        elseif cnt==(3*tN-1)    %last data point of 3rd stg
            intgd_y=zeros(1,kNy);  % data points in dky integral
            for cny=1:kNy   % to calculate integrand for dky integral 
                y=kyRange(cny);
                intgd_x=zeros(1,kNx);  % data points in dkx integral
                for cnx=1:kNx    % to calculate integrand of dkx integral   
                    x=kxRange(cnx);

                    z0=z_fun(x,y,(cnt+1)*dt);
                    dxz=( z_fun(x+derk/2,y,(cnt+1)*dt)-z_fun(x-derk/2,y,(cnt+1)*dt) )/derk;
                    dyz=( z_fun(x,y+derk/2,(cnt+1)*dt)-z_fun(x,y-derk/2,(cnt+1)*dt) )/derk;
                    dtz=( z_fun(x,y,(cnt+1)*dt+dert/2)-z_fun(x,y,(cnt+1)*dt-dert/2) )/dert;
                    if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end

                    intgd_x(cnx)=  2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));                 
                end

                intgd_y(cny)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral    
            end
            discont_t(3)=sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy;   % here I took dky integral, for data point at t=3*T/4 by using H3
        end
    end
    
    % taking time integrals
    T1=sum( intgd_t(0*tN+1:1*tN) + [intgd_t(0*tN+2:1*tN) discont_t(1)] )/2*dt;   % integral within the first time stage, taking care of discontinuity
    T2=sum( intgd_t(1*tN+1:2*tN) + [intgd_t(1*tN+2:2*tN) discont_t(2)] )/2*dt;   % integral within the 2nd time stage..
    T3=sum( intgd_t(2*tN+1:3*tN) + [intgd_t(2*tN+2:3*tN) discont_t(3)] )/2*dt;   % integral within the 3rd time stage..

    intgrl1=T1+T2+T3;
    HopfInvariant1=-1/4/pi^2*intgrl1;
    
%     save(filename,'HopfInvariant1')
% end
toc

HopfInvariant1












branchcut=0;   % 0=>Hopf of zone edge
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------- Return Map -------------------------------------

tic
% if N==1   %no need. this code consider only a single period

intgd_y=zeros(1,kNy);  % data points in dky integral
for cky=1:kNy    % to calculate integrand for dky integral  
    
    intgd_x=zeros(1,kNx);  % data points in dkx integral
    for ckx=1:kNx    % to calculate integrand of dkx integral 
        kx=kxRange(ckx);
        ky=kyRange(cky);
        
                
        zT=U3(kx,ky,T/3)* U2(kx,ky,T/3)* U1(T/3) *z_i;    %states at the end of normal evolution time T
        zx_minus=U3(kx-derk/2,ky,T/3)* U2(kx-derk/2,ky,T/3)* U1(T/3) *z_i; 
        zx_plus=U3(kx+derk/2,ky,T/3)* U2(kx+derk/2,ky,T/3)* U1(T/3) *z_i; 
        zy_minus=U3(kx,ky-derk/2,T/3)* U2(kx,ky-derk/2,T/3)* U1(T/3) *z_i;
        zy_plus=U3(kx,ky+derk/2,T/3)* U2(kx,ky+derk/2,T/3)* U1(T/3) *z_i;
        
        
        [hx,hy,hz, hm]=Heff(kx,ky,branchcut);
        Uret=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        
        [hx,hy,hz, hm]=Heff(kx-derk/2,ky,branchcut);
        Ux_minus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        [hx,hy,hz, hm]=Heff(kx+derk/2,ky,branchcut);
        Ux_plus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        
        [hx,hy,hz, hm]=Heff(kx,ky-derk/2,branchcut);
        Uy_minus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        [hx,hy,hz, hm]=Heff(kx,ky+derk/2,branchcut);
        Uy_plus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);   % inverse time is implemented here !!!!!!!! -T
        
        
        
        intgd_t=zeros(1,3*tN+1);
        %calculate first point t=T here
        dtz=( Uret(dert)*zT - zT )/dert;
        dxz=(zx_plus-zx_minus)/derk;  
        dyz=(zy_plus-zy_minus)/derk;
        intgd_t(1)=2i*zT' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));      
               
        for cnt=1:3*tN 
            z0=Uret(cnt*dt)*zT;
            dtz=( Uret(cnt*dt+dert/2)*zT - Uret(cnt*dt-dert/2)*zT )/dert;            
            dxz=( Ux_plus(cnt*dt)*zx_plus - Ux_minus(cnt*dt)*zx_minus )/derk;
            dyz=( Uy_plus(cnt*dt)*zy_plus - Uy_minus(cnt*dt)*zy_minus )/derk;
            
            if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end
            
            intgd_t(cnt+1)=2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));     
        end
        
        intgd_x(ckx)= sum( intgd_t(1:3*tN)+intgd_t(2:3*tN+1) )/2*dt;
    end
    
    intgd_y(cky)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral   
end

HopfInvariant_ReturnMap=-1/4/pi^2* sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy,


HopfInvariant_pi=HopfInvariant1+HopfInvariant_ReturnMap

% else,  error('code assumes 1period to calculate the return map') 
% end








branchcut=1;   % 1=>Hopf of zone center
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ----------- Return Map -------------------------------------

intgd_y=zeros(1,kNy);  % data points in dky integral
for cky=1:kNy    % to calculate integrand for dky integral  
    
    intgd_x=zeros(1,kNx);  % data points in dkx integral
    for ckx=1:kNx    % to calculate integrand of dkx integral 
        kx=kxRange(ckx);
        ky=kyRange(cky);
        
                
        zT=U3(kx,ky,T/3)* U2(kx,ky,T/3)* U1(T/3) *z_i;    %states at the end of normal evolution time T
        zx_minus=U3(kx-derk/2,ky,T/3)* U2(kx-derk/2,ky,T/3)* U1(T/3) *z_i; 
        zx_plus=U3(kx+derk/2,ky,T/3)* U2(kx+derk/2,ky,T/3)* U1(T/3) *z_i; 
        zy_minus=U3(kx,ky-derk/2,T/3)* U2(kx,ky-derk/2,T/3)* U1(T/3) *z_i;
        zy_plus=U3(kx,ky+derk/2,T/3)* U2(kx,ky+derk/2,T/3)* U1(T/3) *z_i;
        
        
        [hx,hy,hz, hm]=Heff(kx,ky,branchcut);
        Uret=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        
        [hx,hy,hz, hm]=Heff(kx-derk/2,ky,branchcut);
        Ux_minus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        [hx,hy,hz, hm]=Heff(kx+derk/2,ky,branchcut);
        Ux_plus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        
        [hx,hy,hz, hm]=Heff(kx,ky-derk/2,branchcut);
        Uy_minus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);    % inverse time is implemented here !!!!!!!! -T
        
        [hx,hy,hz, hm]=Heff(kx,ky+derk/2,branchcut);
        Uy_plus=@(tt)  (cos(-hm*tt)*eye(2)-1i*sin(-hm*tt)* [hz, hx-1i*hy; hx+1i*hy, -hz]);   % inverse time is implemented here !!!!!!!! -T
        
        
        
        intgd_t=zeros(1,3*tN+1);
        %calculate first point t=T here
        dtz=( Uret(dert)*zT - zT )/dert;
        dxz=(zx_plus-zx_minus)/derk;  
        dyz=(zy_plus-zy_minus)/derk;
        intgd_t(1)=2i*zT' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));      
               
        for cnt=1:3*tN 
            z0=Uret(cnt*dt)*zT;
            dtz=( Uret(cnt*dt+dert/2)*zT - Uret(cnt*dt-dert/2)*zT )/dert;            
            dxz=( Ux_plus(cnt*dt)*zx_plus - Ux_minus(cnt*dt)*zx_minus )/derk;
            dyz=( Uy_plus(cnt*dt)*zy_plus - Uy_minus(cnt*dt)*zy_minus )/derk;
            
            if any([size(z0)~=[2 1] size(dxz)~=[2 1] size(dyz)~=[2 1] size(dtz)~=[2 1]] ),  error('these states are not 2x1 vectors');  end
            
            intgd_t(cnt+1)=2i*z0' * (dxz*imag(dyz'*dtz) + dyz*imag(dtz'*dxz) + dtz*imag(dxz'*dyz));     
        end
        
        intgd_x(ckx)= sum( intgd_t(1:3*tN)+intgd_t(2:3*tN+1) )/2*dt;
    end
    
    intgd_y(cky)=sum(intgd_x(1:kNx-1)+ intgd_x(2:kNx) )/2*dx;   % here I took dkx integral   
end

HopfInvariant_ReturnMap=-1/4/pi^2* sum(intgd_y(1:kNy-1)+ intgd_y(2:kNy))/2*dy,


HopfInvariant_0=HopfInvariant1+HopfInvariant_ReturnMap




toc












% [~,hostname]=system('hostname');
% sendgmail('fnurunal@gmail.com','MatLab:',[ hostname ' found H=' num2str(HopfInvariant) ' over ' num2str(N) 'T. \omega=' num2str(w) ',\delta=' num2str(delta)])
% 












 
function [hx_unit,hy_unit,hz_unit,h_mag] = Heff(kx,ky,branchcut)
global J1 J2 J3 T deltaA 




hz=deltaA;    %fixed at all times
           H1=[hz -J1; -J1 -hz];
           H2=[hz -J2*exp(1i*(kx+sqrt(3)*ky)/2) ; -J2*exp(-1i*(kx+sqrt(3)*ky)/2) -hz];
           H3=[hz -J3*exp(1i*(-kx+sqrt(3)*ky)/2) ; -J3*exp(-1i*(-kx+sqrt(3)*ky)/2) -hz];


          

           %% ----------- Now form Floquet Hamiltonian --------  
           [V1,E1]=eig(H1);   E1=diag(E1);
           [V2,E2]=eig(H2);   E2=diag(E2);
           [V3,E3]=eig(H3);   E3=diag(E3);


           % ----- the basis I'll use is flatbands: [1;0], [0;1] -----------------
           t=T/3;    % length of evolution for each time stage, 4 or 5 stages
           U=zeros(2);    %This will be unitary evolution matrix for time T; at these kx,ky
           Psi1=[1;0];
           Psi2=[0;1];
           % --- now evolve Psi1 ----------------------------------      
           % U1---
           M=V1'*Psi1;    % overlap of Psi with e.st.s of H1, 2-by-1 vector
           C=exp(-1i*E1*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H1                        
           % U2---
           M=(V2'*V1) * C;    % (overlap of e.st.s of H1 with e.st.s of H2) * [c1; c2] = M, 2-by-1
           D=exp(-1i*E2*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H1
           C=D;    %these coefficients are my new C now        
           % U3---
           M=(V3'*V2) * C;    % (overlap of e.st.s of H2 with e.st.s of H3) * [c1; c2] = M, 2-by-1
           D=exp(-1i*E3*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H3
           C=D;    %these coefficients are my new C now    

           A=V3'*Psi1;    % Psi1 in the basis of e.st.s of H3
           B=V3'*Psi2;    % Psi2 in the basis of e.st.s of H3    

           U(1,1)=A'*C;
           U(2,1)=B'*C;


           % --- now evolve Psi2 ----------------------------------      
           % U1---
           M=V1'*Psi2;    % overlap of Psi with e.st.s of H1, 2-by-1 vector
           C=exp(-1i*E1*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H1                        
           % U2---
           M=(V2'*V1) * C;    % (overlap of e.st.s of H1 with e.st.s of H2) * [c1; c2] = M, 2-by-1
           D=exp(-1i*E2*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H1
           C=D;    %these coefficients are my new C now        
           % U3---
           M=(V3'*V2) * C;    % (overlap of e.st.s of H2 with e.st.s of H3) * [c1; c2] = M, 2-by-1
           D=exp(-1i*E3*t) .* M ;    % coefficients, for the evolved state, projected onto e.st.s of H3
           C=D;    %these coefficients are my new C now 

           U(1,2)=A'*C;
           U(2,2)=B'*C;



           [Evec,Eval]=eig(U);   %Diagonalize the time-evolution operator over a period
           En_temp=1i/T*(log(diag(Eval)));   % 2 eigenvalues  

           if max(abs(imag( En_temp)))>1e-7,  putvar(En_temp);   error('this is imaginary');  end
           [~,ind]=max(real(En_temp));   En_temp(ind)=En_temp(ind)-branchcut*2*pi/T; 

           % ------- now construct Heff @kx,ky--------------
           HF= real(En_temp(1))* Evec(:,1)*(Evec(:,1)') + real(En_temp(2))* Evec(:,2)*(Evec(:,2)');

           if max(abs(imag( diag(HF) )))>1e-8, putvar(Evec,HF);  error('diagonals of Heff are imaginary!!!!! and I was just about to delete them!');   end
           HF(1,1)=real(HF(1,1));   HF(2,2)=real(HF(2,2));
           if max(max(HF-HF'))>1e-8,  putvar(Evec,HF);  error('AAAAAAA!!!! Heff not Hermitian. Stillll??!!!');   end

           
           hx=real(HF(2,1));     hy=imag(HF(2,1));     hz=(HF(1,1)-HF(2,2))/2;     h_mag=norm([hx hy hz]);
           hx_unit=hx/h_mag;      hy_unit=hy/h_mag;      hz_unit=hz/h_mag;



















