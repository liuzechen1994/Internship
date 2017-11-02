clc;clear all; close all;
x_b = 0.1; x_e = 0.9; N = 400;
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;

n = 7;
degree = 7;

x_mid = x(101:300);
F = ones(length(x_mid), degree);
scale = (x_mid(end)-x_mid(1))/(2*pi);
mean = (x_mid(1) + x_mid(end))/2;

for i = 2:2:degree
    F(:,i) = cos(1/scale * (i/2) * (x_mid - mean));
end

for i = 3:2:degree
    F(:,i) = sin(1/scale * fix(i/2) * (x_mid - mean));
end

figure()
plot(x_mid,F)
legend('1','2','3','4','5','6','7')

x_front = x(1:100);
x_end = x(301:400);
x_front_end = [x_front x_end];

for i = 0:(n-1)
   B(:,i+1) = nchoosek(n,i).*x_front_end.^i.*(1-x_front_end).^(n-i);
end

% figure()
% plot(x_front_end, B)

F_full = [B(1:100,:); F; B(101:200,:)];
plot(x, F_full)

% A = F_full;
% w = ones(size(A,1),1);
% % w(1:100) = 1e5;
% % w(300:400) = 1e5;
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
% figure()
% plot(x,Q)

for m = 1:size(A,2)
    for n = 1:size(A,2)
        M(m,n) = Q(:,m)'*Q(:,n);
    end
end