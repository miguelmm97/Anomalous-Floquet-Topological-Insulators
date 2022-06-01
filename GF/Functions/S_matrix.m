function [q]=S_matrix(omega,kx,ky,delta,n_t,dt,h_0,J_1,J_2,J_3)
% Calculation of the topological charge with the S-matrix method (PRL Nur
% on how to directly measure Floquete invariants). For a given frecuency (omega),
% point in BZ (kx,ky), potential offset (delta), time resolution (n_t, dt)
% and the offset for the hamiltonian written in the Pauli decomposition (h_0)

change=0.05; % Differentiation infinitesimal
change_t=0.01;

% Differentiation with respect to lambda (omega in this case)
T1 = (2*pi/omega)+change_t; 
t1 = [linspace(0,T1/3,n_t),linspace(T1/3,2*T1/3,n_t),linspace(2*T1/3,T1,n_t)]; % Time vector
dt1 = t1(2)-t1(1); % Resolution of the time vector

T2=(2*pi/omega)-change_t;
t2 = [linspace(0,T2/3,n_t),linspace(T2/3,2*T2/3,n_t),linspace(2*T2/3,T2,n_t)];
dt2 = t2(2)-t2(1);
          
deriv=@(A,B,h) (A-B)/(2*h); % Numerical derivative
                
        

% Spectrum in kx+change
HF=1i*logm(Time_evolution_bulk(delta, kx+change, ky, n_t, dt, J_1, J_2, J_3)*expm(1i*eye(2).*h_0)); % We take the hamiltonian H*T
hx_kx1=real(HF(2,1)); % Pauli component x
hy_kx1=imag(HF(2,1));  % Pauli component y
hz_kx1=(HF(1,1)-HF(2,2))/2; % Pauli component z
% Spectrum in kx-change
HF= 1i*logm(Time_evolution_bulk(delta, kx-change, ky, n_t, dt,J_1,J_2,J_3)*expm(1i*eye(2).*h_0)); %We take the hamiltonian H*T
hx_kx2=real(HF(2,1));% Pauli component x
hy_kx2=imag(HF(2,1));% Pauli component y
hz_kx2=(HF(1,1)-HF(2,2))/2;% Pauli component z


% Spectrum in ky+change
HF= 1i*logm(Time_evolution_bulk(delta, kx, ky+change, n_t, dt,J_1,J_2,J_3)*expm(1i*eye(2).*h_0));
hx_ky1=real(HF(2,1));% Pauli component x
hy_ky1=imag(HF(2,1));% Pauli component y
hz_ky1=(HF(1,1)-HF(2,2))/2;% Pauli component z
% Spectrum in ky-change
HF= 1i*logm(Time_evolution_bulk(delta, kx, ky-change, n_t, dt,J_1,J_2,J_3)*expm(1i*eye(2).*h_0));
hx_ky2=real(HF(2,1));% Pauli component x
hy_ky2=imag(HF(2,1));% Pauli component y
hz_ky2=(HF(1,1)-HF(2,2))/2;% Pauli component z


% Spectrum in T+change
HF= 1i*logm(Time_evolution_bulk(delta, kx, ky, n_t, dt1,J_1,J_2,J_3)*expm(1i*eye(2).*h_0));
hx_T1=real(HF(2,1));% Pauli component x
hy_T1=imag(HF(2,1));% Pauli component y
hz_T1=(HF(1,1)-HF(2,2))/2;% Pauli component z
% Spectrum in T-change
HF= 1i*logm(Time_evolution_bulk(delta, kx, ky, n_t, dt2,J_1,J_2,J_3)*expm(1i*eye(2).*h_0));
hx_T2=real(HF(2,1));% Pauli component x
hy_T2=imag(HF(2,1));% Pauli component y
hz_T2=(HF(1,1)-HF(2,2))/2;% Pauli component z




% Derivatives
d_xx=deriv(hx_kx1, hx_kx2, change);
d_xy=deriv(hx_ky1, hx_ky2, change);
d_xT=deriv(hx_T1, hx_T2, change_t);
d_yx=deriv(hy_kx1, hy_kx2, change);
d_yy=deriv(hy_ky1, hy_ky2, change);
d_yT=deriv(hy_T1, hy_T2, change_t);
d_zx=deriv(hz_kx1, hz_kx2, change);
d_zy=deriv(hz_ky1, hz_ky2, change);
d_zT=deriv(hz_T1, hz_T2, change_t);

% S matrix
S=[d_xx, d_xy, d_xT;
   d_yx, d_yy, d_yT;
   d_zx, d_zy, d_zT];

% Topogical Charge
q=sign(real(det(S)));
  
end
                
                
                