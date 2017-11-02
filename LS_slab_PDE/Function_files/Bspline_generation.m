%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; close all; clear all;
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
% % spatial grid to be used
% N = 400; x_e  = x_meas(end); x_b = x_meas(1);
% dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
% x = x_b:dx:x_e;
% 
% degree  = 4; % d_x
% % choose control points B-spline
% c = linspace(x_b,x_e,8); % n_x
% 
% % Generate intial guesss
% Dc = interp1(xsim,Dsim,c,'pchip'); % ordinates of control points
% Vc = interp1(xsim,Vsim,c,'pchip');
% Kc = interp1(xsim,Ksim,c,'pchip');
% Pc = interp1(xsim,Psim,c,'pchip');
% 
% theta = Dc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [BP,gamma_i,BPc] = Bspline_generation(c,theta,x,degree)
% function [BP,gamma,BPc] = Bspline_generation(c,y,x,degree)
% generates the b-spline matrix and calculates the initial
% fit to the initial guess of the transport
% coefficient  using least-squares fit: gamma = BP\theta
% 
% Input
%   c       control point locations uniformly spaced (C)
%   theta   value of transport coefficient at c, ordinates of control
%           points
%   x       discretized locations (N)
%   degree  B-spline degree (cubic is degree = 4)
%
% Output
%   BP      matrix with B-spline values on x-grid
%   gamma_i B-spline coefficients, i.e., y = BP*gamma
%   BPc     matrix with B-spline values only on c-grid

nr_extra_knots = degree+1; % default
delta_right = (c(end)-c(end-1))*1;
knots_extra_right = linspace(c(end)+delta_right,c(end)+...
                           nr_extra_knots*delta_right,nr_extra_knots);
                       
knots=[c knots_extra_right]; % knots extra right don't understand ask Peter
    
deltac = diff(knots); % Calculate the difference between each knot
knots  = knots-(degree+1)/2.*[deltac deltac(end)];
%% Generate sparse grid B-spline to calculate initial gamma
BPc = 0;
j=0;    
for s_point = c;
    j=j+1;
    for i=length(theta):-1:1,
        BPc(j,i) = BBspline(degree,i,s_point,knots);    
    end       
end,
    
gamma_i = BPc\theta(:);

%% Generate B-splines on finite difference grid for evaluation

BP = 0;
j=0;    
for x_point = x
    j=j+1;
    for i = length(theta):-1:1,
        BP(j,i) = BBspline(degree,i,x_point,knots);
    end       
end

%% Non-uniform rational B-spline
% w = ones(length(c),1);
% % w(1) = 1; w(8) = 0.1;
% w(2:7) = 10;
% 
% M = zeros(length(BP),1);
% for j = 1:length(c)
%     M = M + BP(:,j) * w(j);
% end
% 
% for i = 1:length(c)
%     V(:,i) =  (BP(:,i) * w(i))./M;
% end
    

%%
% A = BP;
% 
% for j  = 1:size(A,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) =Q(:,i)'*A(:,j);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end
% Q(:,1) = 0;

%%
% w = ones(size(A,1),1);
% w(1:20) = 1e2;
% w(390:400) = 1e2;
% 
% 
% for j  = 1:size(A,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) = sum(v.*Q(:,i).*w)/sum(Q(:,i).*Q(:,i).*w);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end
% % Q(:,j) = 0;


%BP = struct('BP',BP);

%% Test
% figure()
% plot(x,BP(:,1:end))
% legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
%       ,'F_8(x)')
% figure()
% plot(V)
% legend('F_1(x)','F_2(x)','F_3(x)','F_4(x)','F_5(x)','F_6(x)','F_7(x)'...
%       ,'F_8(x)')
%% For testing only
% %%
% clear all; clc;
% x = 0:1e-3:2.2; 
% s = linspace(0,2.2,10);
% y = -s.^3+s+5;
% degree = 4;
% 
% y1=BP*gamma; 
% 
% figure
% plot(x,BP,'r',s,BPc,'ko')
% 
% figure
% plot(x,y1,'k',s,BPc*gamma,'o-')
