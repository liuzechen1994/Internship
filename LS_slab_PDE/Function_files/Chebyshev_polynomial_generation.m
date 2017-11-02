function [ T ] = Chebyshev_polynomial_generation( x,degree )
%CHEBYSHEV_POLYNOMIAL_GENERATION Summary of this function goes here
%   Detailed explanation goes here

T  = zeros(length(x),degree);

% scale x to [-1,1]
x_scale = 2/(x(end) - x(1)) * (x - 0.5);
% x_scale = x;

T(:,1) = 1;
T(:,2) = x_scale;

for i = 3:degree
    T(:,i) = 2 * x_scale' .*T(:,i-1) - T(:,i-2);
end

% Normalization
for i = 1:degree
    T(:,i) = T(:,i)/norm(T(:,i));
end

