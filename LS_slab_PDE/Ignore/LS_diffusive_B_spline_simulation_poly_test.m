clear all;
close all;
clc
% %addpath('\\Rijnh\Shares\Users\vanberkel\My Documents\MATLAB\Paper-systematic-optimization_18-04-2016\function_files\\Transfer_function_num')
addpath('\\Rijnh\Shares\Users\vanberkel\My Documents\MATLAB\LS_Bspline_v2\Test_Jacobian')
addpath('\\Rijnh\Shares\Users\vanberkel\My Documents\MATLAB\LS_Bspline_v2\Function_files')
addpath('\\Rijnh\Shares\Users\vanberkel\My Documents\MATLAB\LS_Bspline_v2\Simulation')

set(0,'defaultaxesfontsize',14);
set(0,'defaultlinelinewidth',1.5);

%% ++++++++++++++++++ Generate measurement data +++++++++++++++++++++++++ %%
%% Simulate temperature measurements

%% Input parameters
omega = 2*pi*25; x_b = 0; x_e = 1;
Nsim = 500; xsim = x_b:(x_e-x_b)/(Nsim-1):x_e; rp = 1;
x_meas = linspace(0.1,0.9,20);

% input power definition
Ptot = 0.7; sigma = 0.1; xdep = 0.4; MW2keVs = 1;%6.24e21/2.1e19; % Factor to convert to keV   
Pdep = MW2keVs*Ptot/(sigma*sqrt(pi))*exp(-(xsim-xdep).^2/sigma.^2);

% Values on grid x
Dsim = 10*xsim.^3-xsim+5;
Vsim = -5*xsim.^3-2*xsim-5;
Ksim = 0*ones(size(Dsim));
Psim = Pdep;

% Simulation parameters
% chi = D; Vsim = interp1(x,V,r); tau = interp1(x,K,r);
% Psim =  interp1(x,Pdep,r(1:end-1));
dt= 1e-4; t = 0:dt:0.04-dt; u = square(2*pi*25*t,70); U = fft(u)/length(u);
f = 1;
omega_all = 2*pi*25*f;

[Gh,Rp,y0,profiles] = SlabFD_v2(Dsim,Vsim,Ksim,omega_all,xsim,x_meas,Psim,x_b,x_e,Nsim);
Pomega = U(f+1);
Theta = Gh.*repmat(Pomega(:),1,size(Gh,2));
Gsim = Gh.';
% [yout,tout] = lsim(SYS2,utt,ttt,y0);


% Define input
Input  = [Pomega(:),0*Pomega(:),0*Pomega(:)];
Output =  Theta(:,:);

%
% Output = Theta+0.0000*(randn(size(Theta))+randn(size(Theta))*1i);
% Input = Pomega(:);


%% ++++++++++++++++++++ Input to the code ++++++++++++++++++++++++++++++ %%

%% Settings of signal

% spatial grid to be used
N = 500; x_e  = 1;%x_meas(end); 
x_b = 0;%x_meas(1);
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;

% measurement points

Unknowns = 'DV';






%% Generate intial guesss
Din = interp1(xsim,Dsim,x,'pchip');
Vin = interp1(xsim,Vsim,x,'pchip');
Kin = interp1(xsim,Ksim,x,'pchip');
Pin = interp1(xsim,Psim,x,'pchip');

Pin = ones(size(Pin));
Din = ones(size(Din));
%% Generate B-spline for basis function BF

degree  = 4;
% choose control points B-spline
c = linspace(x_b,x_e,20);

% Generate intial guesss
Dc = interp1(xsim,Dsim,c,'pchip');
Vc = interp1(xsim,Vsim,c,'pchip');
Kc = interp1(xsim,Ksim,c,'pchip');
Pc = interp1(xsim,Psim,c,'pchip');

[inp_mat.BP.D,inp_mat.gamma_i.D,~] = Bspline_generation(c,Dc,x,degree); % Diffusion coefficient
[inp_mat.BP.V,inp_mat.gamma_i.V,~] = Bspline_generation(c,Vc,x,degree); % Convective velocity
[inp_mat.BP.K,inp_mat.gamma_i.K,~] = Bspline_generation(c,Kc,x,degree); % Damping coefficient K_inv
[inp_mat.BP.P,inp_mat.gamma_i.P,~] = Bspline_generation(c,Pc,x,degree); % Power deposition profile

%% Generate polynomial approximation

order = 4;
[inp_mat.gamma_i.D,inp_mat.BP.D] = Polynomial_generation(x,Din,order);
[inp_mat.gamma_i.V,inp_mat.BP.V] = Polynomial_generation(x,Vin,order);
%[inp_mat.gamma_i.P,inp_mat.BP.P] = Polynomial_generation(x,Pin,7);


%% ++++++++++++++++++ Finite difference++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



% Store original profiles
inp_mat.coeff.D = Din; inp_mat.coeff.V = Vin; inp_mat.coeff.K = Kin; inp_mat.coeff.P = Pin;

%% Boundary conditions
% Choose boundary conditions 'Temperature', 'Neumann', 'Diriclet', 'Robin' 


% For now only 'Temperature','Temperature'

% Boundary conditions are measurements
nbc = 0;
xbc_meas = x_meas;
%inp_mat.bc.xa = x_meas(1); inp_mat.bc.xb =  x_meas(end);
xdom = [x_meas(1),xbc_meas,x_meas(end)];

bcl = 1; bcr = 1;
v_bc = 1:N; % determines size of A, B, C matrices
lv_bc = length(v_bc);
%% Sensor locations on grid
C = C_matrix_generation(x,xbc_meas,N); % this is questionable
inp_mat.C = C(:,v_bc);

%% Make complex vector
inp_mat.di = spdiags(1i*ones(lv_bc,1),0,lv_bc,lv_bc);

%% Define Laplacian matrix in terms of D, V, and K_inv
inp_mat.L = Lmatrix_generation(lv_bc,dx);
inp_mat.L.D(1,2) = -inp_mat.L.D(1,1);
inp_mat.L.V(1,2) = -inp_mat.L.V(1,1);


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
%inp_mat.P = Pdep;
inp_mat.v_bc = v_bc;

%%

inp_mat_ideal = inp_mat; inp_mat_ideal.cases.D_on = 0; inp_mat_ideal.cases.V_on = 0;
inp_mat_ideal.cases.K_on = 0; inp_mat_ideal.cases.P_on = 0; 
%inp_mat_ideal.coeff.D = Dest;
gamma_ideal.D = []; gamma_ideal.P = [];

%% Compare transfer functions
[Gideal,~] = Bspline_G_and_J_G_case(inp_mat_ideal,gamma_ideal,omega_all(1));

figure(110)
subplot(2,1,1)
semilogy(x_meas,abs(Gsim(:,1)),'k-',...
         xbc_meas,abs(Gideal(:,1)),'bo')
ylabel('Transfer function')
legend('G_{simulation}','G_{estimation model}','location','NorthWest')
xlabel('x')

grid on
subplot(2,1,2)
semilogy(x_meas,abs(Gsim(:,1)-Gideal(:,1)),'ko-')
ylabel('absolute error')
grid on
xlabel('x')
%% Least-squares scheme
options.maxNbfOfSteps = 300;

[gamma_best, CostBest, CovP] = LM_optimization_Bspline(inp_mat,gamma,omega_all,Output,Input,xdom,options);
CostBest


inp_mat_ideal = inp_mat; inp_mat_ideal.cases.D_on = 0; inp_mat_ideal.cases.V_on = 0;
inp_mat_ideal.cases.K_on = 0; inp_mat_ideal.cases.P_on = 0; 
%inp_mat_ideal.coeff.D = Dest;
gamma_ideal.D = []; gamma_ideal.P = [];
[e] = Cost_DVKP_Bspline_ideal(inp_mat_ideal,gamma_ideal,omega_all,Output,Input,xdom);
CostIdeal = e'*e

%% Compare transfer functions
[Gideal,~] = Bspline_G_and_J_G_case(inp_mat_ideal,gamma_ideal,omega_all(1));
[Gbest,~] = Bspline_G_and_J_G_case(inp_mat,gamma_best,omega_all(1));

figure(111)
subplot(2,1,1)
semilogy(x_meas,abs(Gsim(:,1)),'kx-',...
         xbc_meas,abs(Gideal(:,1)),'bo-',...
         xbc_meas,abs(Gbest(:,1)),'ro-')
grid on
ylabel('Transfer function')
legend('G_{simulation}','G_{estimation model}','G_{best}','location','NorthWest')
xlabel('x')

subplot(2,1,2)
semilogy(xbc_meas,abs((Gideal(:,1)-Gsim(:,1))./Gsim(:,1)),'ko-',...
         xbc_meas,abs((Gideal(:,1)-Gbest(:,1))./Gideal(:,1)),'ro-')
grid on
ylabel('Relative error')
legend('G_{simulation}-G_{estimation model}','G_{best}-G_{estimation model}','location','SouthWest')
xlabel('x')

%%

% plot(xsim,Dsim,'k-o',...
%      x,inp_mat.BP.D*gamma.D,'k--',...
%      x,inp_mat.BP.D*gamma_best.D,'r-')



%%
if inp_mat.cases.D_on
figure
% subplot(2,2,1)
Dmeas = interp1(xsim,Dsim,x_meas);

Dest = inp_mat.BP.D*gamma_best.D;
hold on
plot(xsim,Dsim,'k-',...
     x,inp_mat.BP.D*gamma.D,'k-.',...
     x,inp_mat.BP.D*gamma_best.D,'r--',...
     x_meas,Dmeas,'ko')
 axis([x_meas(1),x_meas(end),0,max(Dsim)])
 box on
 hold off
 grid on
 xlabel('x')
 ylabel('D(x)')
 legend('D simulation','D intial','D estimated','meas. points')
end




if inp_mat.cases.V_on

Vest = inp_mat.BP.V*gamma_best.V;
figure
hold on
plot(xsim,Vsim,'k-',...
     x,inp_mat.BP.V*gamma.V,'k-.',...
     x,inp_mat.BP.V*gamma_best.V,'r--')
 axis([x_meas(1),x_meas(end),min(Vsim),max(Vsim)])
 box on
 stem(x_meas,max(Vest)*ones(size(x_meas)))
 hold off
 grid on
end

if inp_mat.cases.K_on
subplot(2,2,3)
plot(x,K,'k-',...
     x,inp_mat.BP.K*gamma.K,'k-.',...
     x,inp_mat.BP.K*gamma_best.K,'r--',...
     c,BPc.K*gamma_best.K,'bo')
end

if inp_mat.cases.P_on
    Pest = inp_mat.BP.P*gamma_best.P;
    Pmeas = interp1(xsim,Pdep,x_meas);

figure
plot(xsim,Pdep,'k-',...
     x,Pin,'k-.',...
     x,Pest,'r--',...
     x_meas,Pmeas,'ko')
grid on


 xlabel('x')
 ylabel('P(x)')
 legend('P simulation','P intial','P estimated','meas. points')
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

