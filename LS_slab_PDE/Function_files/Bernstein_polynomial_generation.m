function [Q,cof ] = Bernstein_polynomial_generation(x,theta,n)
%BERNSTEIN_POLYNOMIAL_GENERATION Summary of this function goes here
%   Detailed explanation goes here

% Generate Bernstein polynomial
for i = 0:n
    B(:,i+1) = nchoosek(n,i).*x.^i.*(1-x).^(n-i);
end

% Othonormalization
A = B;

w = ones(size(B,1),1);
% w(1:100) = 1e5;
% w(300:400) = 1e5;
 
for j  = 1:size(B,2)
    v = A(:,j);
    for i = 1:j-1
        R(i,j) = sum(v.*Q(:,i).*w)/sum(Q(:,i).*Q(:,i).*w);
        v = v - R(i,j)*Q(:,i);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

cof = Q\(theta)';

end

