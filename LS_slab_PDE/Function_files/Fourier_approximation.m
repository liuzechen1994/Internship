function [Q,cof ] = Fourier_approximation(x,theta,degree)
%FOURIER_APPROXIMATION Summary of this function goes here
%   Detailed explanation goes here

F = ones(length(x), degree);
scale = (x(end)-x(1))/(2*pi);

for i = 2:2:degree
    F(:,i) = cos(1/scale * (i/2) * (x - 0.5));
end

for i = 3:2:degree
    F(:,i) = sin(1/scale * fix(i/2) *(x - 0.5));
end

% theta = theta';
%  cof = (theta\F)';
%  cof = 0.5*ones(7,1)+3;
 
%% Initial 
A = F;
%% Classical Gram-Schmidt 

% for j  = 1:size(A,2)
%     v = A(:,j);
%     for i = 1:j-1
%         R(i,j) =Q(:,i)'*A(:,j);
%         v = v - R(i,j)*Q(:,i);
%     end
%     R(j,j) = norm(v);
%     Q(:,j) = v/R(j,j);
% end

%% Weight Gram-Schmidt orthonormalization

 w = ones(size(F,1),1);
 w(1:20) = 1e2;
 w(390:400) = 1e2;


for j  = 1:size(F,2)
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = sum(v.*Q(:,i).*w)/sum(Q(:,i).*Q(:,i).*w);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end
 
% Q(:,1) = 0;
% Q(:,1) = Q(:,1).*(1/10000);
% Q(:,1) = Q(:,1)/norm(Q(:,1));


 cof = Q\(theta)' + 5;

end

