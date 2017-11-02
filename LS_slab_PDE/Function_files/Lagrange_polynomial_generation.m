function [Q] = Lagrange_polynomial_generation(c,theta,x,num_p)
% LAGRANGE_POLYNOMIAL_GENERATION generate basis function and gamma
% Input
%   c      control points
%   theta  value of transport coefficient
%   x      discretized locations
%   num_p  number of points 

% Output
%   l      matrix with lagrange polynomial on x-grid
%   gamma  cofficients


l = ones(num_p, length(x));
for i = 1:num_p
    for j = 1:num_p
        if (i~=j)
            l(i,:) = l(i,:).* (x - c(j))/(c(i)-c(j));
        end
    end
end
l = l';

% gamma = (theta\l)';
% Gram-Schmidt orthonormalization
A = l;

for j  = 1:size(A,2)
    v = A(:,j);
    for i = 1:j-1
        R(i,j) =Q(:,i)'*A(:,j);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end
%%
% w = ones(size(A,1),1);
%  w(1:20) = 10;
%  w(390:400) = 10;
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
% gamma = Q\(theta)';
end

