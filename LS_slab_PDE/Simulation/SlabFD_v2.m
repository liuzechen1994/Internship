%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parameters
% omega = 2*pi*25; x_b = 0; x_e = 1;
% Nsim = 400; xsim = x_b:(x_e-x_b)/(Nsim-1):x_e; rp = 1;
% x_meas = linspace(0.1,0.9,10);
% 
% % input power definition
% Ptot = 0.7; sigma = 0.1; xdep = 0.5; MW2keVs = 1;%6.24e21/2.1e19; % Factor to convert to keV   
% Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep).^2/sigma.^2); % Gaussian Function
% 
% 
% % Values on grid x
% % First estimate P(x) and D(x)
% Dsim = 10*xsim.^3-xsim+5;
% Vsim = 0*xsim;
% Ksim = 0*ones(size(Dsim));
% Psim = Pdep;
% 
% % Simulation parameters
% % chi = D; Vsim = interp1(x,V,r); tau = interp1(x,K,r); 
% % Psim =  interp1(x,Pdep,r(1:end-1));
% dt= 1e-4; t = 0:dt:0.04-dt; u = square(2*pi*25*t,70); U = fft(u)/length(u); % u = p(t); U = p(omega)
% f = 1:4;
% omega_all = 2*pi*25*f;
% 
% D_x = Dsim;
% V_x = Vsim;
% K_x = Ksim;
% 
% N = Nsim;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ghnum,Rp,y0,profiles] = SlabFD_v2(D_x,V_x,K_x,omega_all,xsim,x_meas,Pdep,x_b,x_e,N)

% This function calculates the transfer functions from P -> Te(X) of the PDE:
  % dT/dt = chi*dT^2/dx^2 + V*dT/dx + tau*T + Pdep
%   X = location of the sensors
%   xe = location of the outer edge
%   N = number of discretization points
%   rp = radial coord. of the profile given
%   rpd = radial coord. of the deposition profile
%   chi = diffusion coefficient/profile
%   V  = convection coefficient/profile
%   tau = damping coefficient/profile

%% +++++++++++++++++++ Generate discretization +++++++++

dx = (x_e-x_b)/(N-1);  % dx = 0.0025
xsim = x_b:dx:x_e;

% Check if profiles are given and interpolate

% if length(rp) == 1;
%     Chi_e = chi+0*x;
%     V_e = V+0*x;
%     tau_inv = tau+0*x;
% elseif norm(rp - x)==0;
%     Chi_e = chi;
%     V_e = V;
%     tau_inv = tau;
% else
%     Chi_e = interp1(rp,chi,x,'cubic');
%     V_e = interp1(rp,V,x,'cubic'); 
%     tau_inv = interp1(rp,tau,x,'cubic');
% end


% interpolate source
% if norm(rpd - x)==0;
%     P = Pdep;
% else P = interp1(rpd,Pdep,x,'cubic');
% end

%% +++++++++++++++++++ A-matrix ++++++++++++++++++++++++
% Generate coefficients    
a_r = D_x/(dx^2) + (V_x)/(2*dx);
b_r = -2*D_x/(dx^2) + K_x; % £¿
c_r = D_x/(dx^2) - (V_x)/(2*dx);

% Generate matrix
%A = spdiags([c_r' b_r' a_r'],-1:1,N,N).'; % ?
A = spdiags([c_r' b_r' a_r'],-1:1,N-1,N-1).'; % ?

% Forward finite difference
A(1,1) = D_x(1)/(dx^2) - V_x(1)/(dx) + K_x(1) ;
A(1,2) = -2*D_x(1)/(dx^2) + V_x(1)/(dx);
A(1,3) = D_x(1)/(dx^2);

% Backward finite difference
% A(N,N-2) = D_x(N)/(dx^2);
% A(N,N-1) = -2*D_x(N)/(dx^2) - V_x(N)/(dx);
% A(N,N) = D_x(N)/(dx^2) + V_x(N)/(dx) + K_x(N);

% Boundary conditions: gpunt= (g(x+h)-g(x-h))/2h=0
% dT(x_1)/dx = 0, T(x_n) = 0
% A(1,2)   = a_r(1)+c_r(1); % ?
A(1,2) = -c_r(1); % ?
A(1,3) = 0;

%% Input matrix
% B(:,1) = Pdep; % deposition profile
% B(N,2) = 0; % Temperature at boundary (not relevant for TF) %
B = Pdep(1:end-1)';
B(N-1,2) = 0; 

%% +++++++++++++++++++ C-matrix ++++++++++++++++++++++++

% Select source and measurement locations
Sensor = interp1(xsim,xsim,x_meas,'nearest');
for i = length(Sensor):-1:1;
    indSf(i) = find(xsim==Sensor(i)); % Find points' number in x that is the same as ith number of Sensor
end
if indSf(end) == length(xsim);
    indSf = indSf(1:end-1); % Boundry condition T_i+1 = 0,no need to put sensor on it
end

% Measurement locations
Rp = xsim(indSf); % Measurement positions

% C = zeros(length(Sensor),N); j = 1;
C = zeros(length(Sensor),N-1); j = 1;
for i = indSf 
    C(j,i) = 1;
    j=1+j;
end
C = sparse(C); %S = sparse(A) converts a full matrix into sparse form by squeezing out any zero elements. If a matrix contains many zeros, converting the matrix to sparse storage saves memory
% D = zeros(length(vector),2);
%% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% +++++++++++++++++++ Calculate transfer function++++++++++++++++++++

%di = spdiags(1i*ones(N,1),0,N,N); %diagonal matrix
di = spdiags(1i*ones(N,1),0,N-1,N-1); %diagonal matrix
for j = length(omega_all):-1:1
    sI_A = omega_all(j)*di-A;
    Ghnum(j,:) = (C/sI_A)*B(:,1);
end

% Steady-state value
y0 = -C/A*B; %0 = x_dot = Ax+Bu

profiles = [xsim(:),D_x(:),V_x(:),K_x(:),Pdep(:)];
%  end

