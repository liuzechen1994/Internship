clear all;
close all;
clc
% %addpath('\\Rijnh\Shares\Users\vanberkel\My Documents\MATLAB\Paper-systematic-optimization_18-04-2016\function_files\\Transfer_function_num')
addpath('C:\Users\ux501\Dropbox\Internship\LS_slab_PDE\Test_Jacobian')
addpath('C:\Users\ux501\Dropbox\Internship\LS_slab_PDE\Function_files')
addpath('C:\Users\ux501\Dropbox\Internship\LS_slab_PDE\Simulation')

set(0,'defaultaxesfontsize',14);
set(0,'defaultlinelinewidth',1.5);

%% ++++++++++++++++++ Generate measurement data +++++++++++++++++++++++++ %%
%% Simulate temperature measurements

% Input parameters
omega = 2*pi*25; x_b = 0; x_e = 1;
Nsim = 400; xsim = x_b:(x_e-x_b)/(Nsim-1):x_e; rp = 1;
x_meas = linspace(0.1,0.9,10);

% n_s = cos(((2*(1:10)-1)/20)*pi);
% x_meas = (n_s + 1)*0.4 + 0.1;

% x_meas(1:4) = linspace(0.1,0.3,4);
% x_meas(5:10) = linspace(0.3,0.7,6);
% x_meas(11:14) = linspace(0.7,0.9,4);

% input power definition
Ptot = 0.7; sigma = 0.1; xdep1 = 0.5; xdep2 = 0.8; MW2keVs = 1;%6.24e21/2.1e19; % Factor to convert to keV   
Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep1).^2/sigma.^2); % Gaussian Function
% Pdep2 = 0.3*MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep2).^2/sigma.^2);
% Pdep = Pdep1 + Pdep2;


% Values on grid x
% First estimate P(x) and D(x)
Dsim = 10*xsim.^3-3*xsim+5;
Vsim = 0*xsim;
Ksim = 0*ones(size(Dsim));
Psim = Pdep;

% Simulation parameters
% chi = D; Vsim = interp1(x,V,r); tau = interp1(x,K,r); 
% Psim =  interp1(x,Pdep,r(1:end-1));
dt= 1e-4; t = 0:dt:0.04-dt; u = square(2*pi*25*t,70); U = fft(u)/length(u); % u = p(t); U = p(omega)
f = 1:4;
omega_all = 2*pi*25*f;
    
[Gh,Rp,y0,profiles] = SlabFD_v2(Dsim,Vsim,Ksim,omega_all,xsim,x_meas,Pdep,x_b,x_e,Nsim); % p8 calculate G_sim
Pomega = U(f+1); % Dc value
Theta = Gh.*repmat(Pomega(:),1,size(Gh,2)); % p9 equation

% [yout,tout] = lsim(SYS2,utt,ttt,y0);


% Define input
Input  = [Pomega(:),Theta(:,1),Theta(:,end)]; % ?
Output =  Theta(:,2:end-1);

%
% Output = Theta+0.0001*(randn(size(Theta))+randn(size(Theta))*1i);
% Input = Pomega(:);


%% ++++++++++++++++++++ Input to the code ++++++++++++++++++++++++++++++ %%

%% Settings of signal

% spatial grid to be used
N = 400; x_e  = x_meas(end); x_b = x_meas(1);
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;

% measurement points

Unknowns = 'P'; 









%% Generate B-spline for basis function BF

degree_p  = 4; % d_x
% choose control points B-spline
c = linspace(x_b,x_e,8); % control points, n_x = 8

% Generate intial guesss
Dc = interp1(xsim,Dsim,c,'pchip'); 
Vc = interp1(xsim,Vsim,c,'pchip');
Kc = interp1(xsim,Ksim,c,'pchip');
Pc = interp1(xsim,Psim,c,'pchip');

%Pc = ones(size(interp1(xsim,Psim,c,'pchip')));

% inp_mat.BP basis function, inp_mat.gamma_i coefficients
[inp_mat.BP.D,inp_mat.gamma_i.D,~] = Bspline_generation(c,Dc,x,degree_p); % Diffusion coefficient
[inp_mat.BP.V,inp_mat.gamma_i.V,~] = Bspline_generation(c,Vc,x,degree_p); % Convective velocity
[inp_mat.BP.K,inp_mat.gamma_i.K,~] = Bspline_generation(c,Kc,x,degree_p); % Damping coefficient K_inv
[inp_mat.BP.P,inp_mat.gamma_i.P,~] = Bspline_generation(c,Pc,x,degree_p); % Power deposition profile

%% Generate polynomial approximation

% Bpoly = [x.^5;x.^4;x.^3;x.^2;x.^1;x.^0]; inp_mat.BP.D = Bpoly.'; inp_mat.gamma_i.D = 0.1*ones(size(inp_mat.BP.D,2),1);
% inp_mat.gamma_i.D = [0;1;10;0;-1;5]; %Radom value or ?

% Bpoly = [x.^1;x.^0]; inp_mat.BP.D = Bpoly.'; inp_mat.gamma_i.D = 10*ones(size(inp_mat.BP.D,2),1);
% inp_mat.gamma_i.D = [0;10];

%% Generate Chebyshev polynomial approximation

degree = 4;
 
% [inp_mat.BP.D] = Chebyshev_polynomial_generation(x,degree);
% inp_mat.gamma_i.D = [150;35;20;5];
% inp_mat.gamma_i.D = [1;0;0;0;0;0];
%% Generate Lagrange polynomial approximation

num_p_D = 4; 
num_p_P = 8;
c_D = linspace(x_b,x_e,num_p_D);
c_P = linspace(x_b,x_e,num_p_P);
Dc = interp1(xsim,Dsim,c_D,'pchip');
Pc = interp1(xsim,Psim,c_P,'pchip');

[inp_mat.BP.D] = Lagrange_polynomial_generation(c_D,Dc,x,num_p_D);
inp_mat.gamma_i.D = [50; 30; 100; 65];
% [inp_mat.BP.P,inp_mat.gamma_i.P] = Lagrange_polynomial_generation(c_P,Pc,x,num_p_P);

%% Generate Fourier approximation

degree_P = 7;
degree_D = 7;
Dc = interp1(xsim,Dsim,x,'pchip');
Pc = interp1(xsim,Psim,x,'pchip');

 [inp_mat.BP.P,inp_mat.gamma_i.P] = Fourier_approximation(x,Pc,degree_P);
[inp_mat.BP.D,inp_mat.gamma_i.D] = Fourier_approximation(x,Dc,degree_D);

%% Generate Bernstein approximation
n = 10;
Pc = interp1(xsim,Psim,x,'pchip');

% [inp_mat.BP.P,inp_mat.gamma_i.P] = Bernstein_polynomial_generation(x,Pc,n);

%% Generate Combined approximation
degree = 7;
Pc = interp1(xsim,Psim,x,'pchip');

% [inp_mat.BP.P,inp_mat.gamma_i.P] = Combined_basis_function(x,Pc,degree);
%% ++++++++++++++++++ Finite difference++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Generate intial guesss
Din = interp1(xsim,Dsim,x,'pchip');
Vin = interp1(xsim,Vsim,x,'pchip');
Kin = interp1(xsim,Ksim,x,'pchip');
Pin = interp1(xsim,Psim,x,'pchip');


% Store original profiles
inp_mat.coeff.D = Din; inp_mat.coeff.V = Vin; inp_mat.coeff.K = Kin; inp_mat.coeff.P = Pin;

%% Boundary conditions
% Choose boundary conditions 'Temperature', 'Neumann', 'Diriclet', 'Robin' 


% For now only 'Temperature','Temperature'

% Boundary conditions are measurements
nbc = 2;
xbc_meas = x_meas(2:end-1);
%inp_mat.bc.xa = x_meas(1); inp_mat.bc.xb =  x_meas(end);
xdom = [x_meas(1),xbc_meas,x_meas(end)];

bcl = 1; bcr = 1;
v_bc = 2:N-1; % determines size of A, B, C matrices
lv_bc = length(v_bc);

%% Sensor locations on grid
C = C_matrix_generation(x,xbc_meas,N); % this is questionable
inp_mat.C = C(:,v_bc);


%% Make complex vector
inp_mat.di = spdiags(1i*ones(lv_bc,1),0,lv_bc,lv_bc);

%% Define Laplacian matrix in terms of D, V, and K_inv
inp_mat.L = Lmatrix_generation(lv_bc,dx);
%inp_mat.L.D(1,2) = -inp_mat.L.D(1,1);

inp_mat.cases = Cases_variables_estimated(Unknowns);
% ==================================================================
% ==================================================================
%% +++++++++++++++++++ Algorithm starts here +++++++++++++++++++
if inp_mat.cases.D_on; gamma.D = inp_mat.gamma_i.D;
else gamma.D = []; end
    
if inp_mat.cases.V_on; gamma.V = inp_mat.gamma_i.V;
else gamma.V = []; end

if inp_mat.cases.K_on; gamma.K = inp_mat.gamma_i.K;
else gamma.K = []; end

if inp_mat.cases.P_on; gamma.P = inp_mat.gamma_i.P;
else gamma.P = []; end


%% Select parameters to be estimated


inp_mat.dx = dx;
inp_mat.P = Pdep;
inp_mat.v_bc = v_bc;




%% Least-squares scheme
options.maxNbfOfSteps = 500;

[gamma_best, CostBest, CovP] = LM_optimization_Bspline(inp_mat,gamma,omega_all,Output,Input,xdom,options);
CostBest
%Dest = inp_mat.BP.D*gamma_best.D;




inp_mat_ideal = inp_mat; inp_mat_ideal.cases.D_on = 0; inp_mat.cases.V_on = 0;
inp_mat_ideal.cases.K_on = 0; inp_mat_ideal.cases.P_on = 0; 
%inp_mat_ideal.coeff.D = Dest;
gamma_ideal.D = []; gamma_ideal.P = [];
[e] = Cost_DVKP_Bspline_ideal(inp_mat_ideal,gamma_ideal,omega_all,Output,Input,xdom);
CostIdeal = e'*e


gamma_initial_P = inp_mat.gamma_i.P
gamma_best_P = gamma_best.P
gamma_ideal_P = inp_mat.BP.P\Pdep'

% gamma_initial_D = inp_mat.gamma_i.D
% gamma_best_D = gamma_best.D
% gamma_ideal_D = inp_mat.BP.D\(Dsim')
%%
if inp_mat.cases.D_on
figure
%subplot(2,2,1)
Dest = inp_mat.BP.D*gamma_best.D;
hold on
plot(xsim,Dsim,'k-',...
     x,inp_mat.BP.D*gamma.D,'b--',...
     x,inp_mat.BP.D*gamma_best.D,'r-d')
plot(x, inp_mat.BP.D*gamma_ideal_D,'g-.')
 legend('D simulation','D initial','D estimated','D ideal')
 xlabel('x')
 ylabel('D(x)')
 %axis([0 1 0 14]);
 box on
 stem(x_meas,max(Dest)*ones(size(x_meas)))
 hold off
 grid on
end
if inp_mat.cases.V_on
figure
plot(rpd,Vsim,'k',...
     x,inp_mat.BP.V*gamma.V,'k.',...
     x,inp_mat.BP.V*gamma_best.V,'r--')
end

if inp_mat.cases.K_on
subplot(2,2,3)
plot(x,K,'k',...
     x,inp_mat.BP.K*gamma.K,'k.',...
     x,inp_mat.BP.K*gamma_best.K,'r--',...
     c,BPc.K*gamma_best.K,'bo')
end
if inp_mat.cases.P_on
    Pest = inp_mat.BP.P*gamma_best.P;
figure
% subplot(2,1,2)
hold on
plot(xsim,Pdep,'k-',...
     x,inp_mat.BP.P*gamma.P,'b--',...
     x,Pest,'r-o')
plot(x, inp_mat.BP.P*gamma_ideal_P,'g-.')
legend('P simulation','P initial','P estimated','P ideal')
xlabel('x')
ylabel('P(x)')
axis([0 1 -2 4]);
grid on
stem(x_meas,max(Pest)*ones(size(x_meas)))
hold off

end

%% Test
figure()
plot(x,inp_mat.BP.P(:,1:end))
grid on
tit = title('Basis function P(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');
xlim([0.1 0.9])
%axis([0.1 0.9 -1 1])
legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
      ,'F_8(x)','F_9(x)','F_{10}(x)','F_{11}(x)')

% figure()
% plot(x,inp_mat.BP.D(:,1:end))
% grid on
% tit = title('Fourier basis function D(x)');
% %set(tit,'interpreter','latex');
% xlabel('x');
% ylabel('y');
% xlim([0.1 0.9])
% %axis([0.1 0.9 -1 1])
% legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
%       ,'F_8(x)','F_9(x)')  
  
% Test orthonormal
for m = 1:size(inp_mat.BP.P,2)
    for n = 1:size(inp_mat.BP.P,2)
        M(m,n) = inp_mat.BP.P(:,m)' * inp_mat.BP.P(:,n);
    end
end

for m = 1:size(inp_mat.BP.D,2)
    for n = 1:size(inp_mat.BP.D,2)
        N(m,n) = inp_mat.BP.D(:,m)' * inp_mat.BP.D(:,n);
    end
end
return

%% ++++++++++++++++++++++++++ Test Jacobian +++++++++++++++++++++++++++ %%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

[Gout,Jout] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);

dtheta = 1e-8; 

   
%% Diffusion
if inp_mat.cases.D_on;
Jnum1D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,1);
Jnum2D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,2);
Jnum3D = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,3);
JD_error = [dtheta,norm(Jnum1D-Jout.Jd(:,:,1)),norm(Jnum2D-Jout.Jd(:,:,2)),norm(Jnum3D-Jout.Jd(:,:,3))]
end
%% Convection
if inp_mat.cases.V_on; 
Jnum1V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,1);
Jnum2V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,2);
Jnum3V = Test_Jacobian_V(inp_mat,gamma,omega,dtheta,3);
JV_error = [dtheta,norm(Jnum1V-Jout.Jv(:,:,1)),norm(Jnum2V-Jout.Jv(:,:,2)),norm(Jnum3V-Jout.Jv(:,:,3))]
end
%% Damping
if inp_mat.cases.K_on; 
Jnum1K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,1);
Jnum2K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,2);
Jnum3K = Test_Jacobian_K(inp_mat,gamma,omega,dtheta,3);
JK_error = [dtheta,norm(Jnum1K-Jout.Jk(:,:,1)),norm(Jnum2K-Jout.Jk(:,:,2)),norm(Jnum3K-Jout.Jk(:,:,3))]
end
%% P
if inp_mat.cases.P_on; 
JnumP = Test_Jacobian_P(inp_mat,gamma,omega,dtheta,1);
JP_error = [dtheta,norm(JnumP-Jout.Jp(:,:,1))]
end

%% % Test full Jacobian
[E,JE] =  Cost_DVKP_Bspline(inp_mat,gamma,omega_all,Output,[Input,0*Input,0*Input],xdom);

[Jnum] = Test_Jacobian_full(inp_mat,gamma,omega_all,Output,[Input,0*Input,0*Input],xdom,dtheta);

Jfull_error = norm(JE-Jnum)


%% Check Orthonormal

for m = 1:size(inp_mat.BP.P,2)
    for n = 1:size(inp_mat.BP.P,2)
        M(m,n) = inp_mat.BP.P(:,m)' * inp_mat.BP.P(:,n);
    end
end
%% Figure presentation
figure(123)
hold on
plot(xsim,inp_mat.BP.P*inp_mat.gamma_i.P,'k.-')
for i = 1:length(gamma.P);
    plot(xsim,inp_mat.BP.P(:,i).*gamma.P(i))
end
plot(xsim,inp_mat.BP.P*inp_mat.gamma_i.P,'k--')
hold off
box on
grid on
xlabel('x')
ylabel('P(x)')
legend('P(x_i)','B-splines')
%% Basis function plot
figure(1)
plot(x,inp_mat.BP.P(:,1:end))
grid on
tit = title('Fourier basis function P(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');
xlim([0.1 0.9])
%axis([0.1 0.9 -1 1])
legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
      ,'F_8(x)','F_9(x)','F_{10}(x)','F_{11}(x)')
  
figure(2)
plot(x,inp_mat.BP.D(:,1:end))
grid on
tit = title('Fourier basis function D(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');
axis([0.1 0.9])
%axis([0.1 0.9 -1 1])
legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
      ,'F_8(x)','F_9(x)')
  
figure(3)
plot(x,inp_mat.BP.D(:,1:end))
grid on
tit = title('Chebyshev basis function D(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');

axis([0.1 0.9 -1 1])
legend('T_1(x)','T_2(x)','T_3(x)','T_4(x)','T_5(x)')

figure(4)
plot(x,inp_mat.BP.D(:,1:end))
grid on
tit = title('Lagrange basis function D(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');

axis([0.1 0.9 -1 1])
legend('L_1(x)','L_2(x)','L_3(x)','L_4(x)')

figure(5)
plot(x,inp_mat.BP.P(:,1:end))
grid on
tit = title('Lagrange basis function P(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');

axis([0.1 0.9 -1 1])
legend('L_1(x)','L_2(x)','L_3(x)','L_4(x)','L_5(x)','L_6(x)','L_7(x)'...
      ,'L_8(x)')

figure(5)
plot(x,inp_mat.BP.P(:,1:end))
grid on
tit = title('B-spline basis function P(x)');
%set(tit,'interpreter','latex');
xlabel('x');
ylabel('y');

axis([0.1 0.9 -1 1])
legend