function [q]=S_matrix_lambda(lambda,kx,ky,delta,n_t,dt,h_0,J_1,J_2,J_3)
% Calculation of the topological charge with the S-matrix method (PRL Nur
% on how to directly measure Floquete invariants). For a given frecuency (omega),
% point in BZ (kx,ky), potential offset (delta), time resolution (n_t, dt)
% and the offset for the hamiltonian written in the Pauli decomposition (h_0)

change=0.05; % Differentiation infinitesimal         
deriv=@(A,B,h) (A-B)/(2*h); % Numerical derivative
                
       
% Spectrum in kx+change
HF=1i*logm(Time_evolution_bulk_lambda(delta, kx+change, ky, n_t, dt, J_1, J_2, J_3, lambda)*expm(1i*eye(2).*h_0)); % We take the hamiltonian H*T
hx_kx1=real(HF(2,1)); % Pauli component x
hy_kx1=imag(HF(2,1));  % Pauli component y
hz_kx1=(HF(1,1)-HF(2,2))/2; % Pauli component z
% Spectrum in kx-change
HF= 1i*logm(Time_evolution_bulk_lambda(delta, kx-change, ky, n_t, dt,J_1,J_2,J_3, lambda)*expm(1i*eye(2).*h_0)); %We take the hamiltonian H*T
hx_kx2=real(HF(2,1));% Pauli component x
hy_kx2=imag(HF(2,1));% Pauli component y
hz_kx2=(HF(1,1)-HF(2,2))/2;% Pauli component z


% Spectrum in ky+change
HF= 1i*logm(Time_evolution_bulk_lambda(delta, kx, ky+change, n_t, dt,J_1,J_2,J_3, lambda)*expm(1i*eye(2).*h_0));
hx_ky1=real(HF(2,1));% Pauli component x
hy_ky1=imag(HF(2,1));% Pauli component y
hz_ky1=(HF(1,1)-HF(2,2))/2;% Pauli component z
% Spectrum in ky-change
HF= 1i*logm(Time_evolution_bulk_lambda(delta, kx, ky-change, n_t, dt,J_1,J_2,J_3, lambda)*expm(1i*eye(2).*h_0));
hx_ky2=real(HF(2,1));% Pauli component x
hy_ky2=imag(HF(2,1));% Pauli component y
hz_ky2=(HF(1,1)-HF(2,2))/2;% Pauli component z


% Spectrum in delta+change
HF= 1i*logm(Time_evolution_bulk_lambda(delta, kx, ky, n_t, dt,J_1,J_2,J_3, lambda+change)*expm(1i*eye(2).*h_0));
hx_d1=real(HF(2,1));% Pauli component x
hy_d1=imag(HF(2,1));% Pauli component y
hz_d1=(HF(1,1)-HF(2,2))/2;% Pauli component z
% Spectrum in delta-change
HF= 1i*logm(Time_evolution_bulk_lambda(delta, kx, ky, n_t, dt,J_1,J_2,J_3, lambda-change)*expm(1i*eye(2).*h_0));
hx_d2=real(HF(2,1));% Pauli component x
hy_d2=imag(HF(2,1));% Pauli component y
hz_d2=(HF(1,1)-HF(2,2))/2;% Pauli component z




% Derivatives
d_xx=deriv(hx_kx1, hx_kx2, change);
d_xy=deriv(hx_ky1, hx_ky2, change);
d_xd=deriv(hx_d1, hx_d2, change);
d_yx=deriv(hy_kx1, hy_kx2, change);
d_yy=deriv(hy_ky1, hy_ky2, change);
d_yd=deriv(hy_d1, hy_d2, change);
d_zx=deriv(hz_kx1, hz_kx2, change);
d_zy=deriv(hz_ky1, hz_ky2, change);
d_zd=deriv(hz_d1, hz_d2, change);

% S matrix
S=[d_xx, d_xy, d_xd;
   d_yx, d_yy, d_yd;
   d_zx, d_zy, d_zd];

% Topogical Charge
q=sign(real(det(S))); 
  
end
                