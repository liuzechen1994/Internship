%%
clc; clear all;
x_b = 0; x_e = 1;
Nsim = 400; xsim = x_b:(x_e-x_b)/(Nsim-1):x_e;
Dsim = 10*xsim.^3-xsim+5;

%%
x_b = 0.1; x_e = 0.9; N = 400;
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;

degree = 7;
Dc = interp1(xsim,Dsim,x,'pchip');

F = ones(length(x), degree);
scale = (x(end)-x(1))/(2*pi);
mean = (x(1) + x(end))/2;

for i = 2:2:degree
    F(:,i) = cos(1/scale * (i/2) * (x - mean));
end

for i = 3:2:degree
    F(:,i) = sin(1/scale * fix(i/2) *(x - mean));
end

%% Gram-Schmidt orthonormalization
A = F;

% Classical Gram-Schmidt
% for j  = 1:size(A,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) =Q(:,i)'*A(:,j);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end

 w = ones(size(F,1),1);
% w(1:100) = 1e5;
% w(300:400) = 1e5;
 
for j  = 1:size(F,2)
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = sum(v.*Q(:,i).*w)/sum(Q(:,i).*Q(:,i).*w);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

%Q(:,1) = 0;

%%
% [m, n] = size(A);
% Q = zeros(m, n, 'like',A);
% R = zeros(n, n, 'like',A);
% QQ = zeros(m, n, 'like',A);
% 
% for j = 1:n
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) = (v.' * QQ(:,i));
%         v = v - R(i,j) * Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v / R(j,j);
%     QQ(:,j) = (conj(Q(:,j)) .* w) ./ (w.' * (Q(:,j).*conj(Q(:,j))));
% end

% 
% 
% 
%% Check orthonormal

for m = 1:size(F,2)
    for n = 1:size(F,2)
        M(m,n) = Q(:,m)' * Q(:,n);
    end
end

figure()
plot(x,F(:,1:end))

figure()
plot(x,Q(:,1:end))

% Dc = Dc';
% cof = (Dc\A)';
