function [ F_full,cof ] = Combined_basis_function( x,theta,degree )
%COMBINED_BASIS_FUNCTION Summary of this function goes here
%   Detailed explanation goes here
n = degree;

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

x_front = x(1:100);
x_end = x(301:400);
x_front_end = [x_front x_end];

for i = 0:(n-1)
   B(:,i+1) = nchoosek(n,i).*x_front_end.^i.*(1-x_front_end).^(n-i);
end

F_full = [B(1:100,:); F; B(101:200,:)];

cof = F_full\(theta)';
end

