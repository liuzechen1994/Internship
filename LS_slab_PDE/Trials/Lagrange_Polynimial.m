x_b = 0; x_e = 1;
Nsim = 400; xsim = x_b:(x_e-x_b)/(Nsim-1):x_e;
Dsim = 10*xsim.^3-xsim+5;
%%
N = 400; x_e  = 0.9; x_b = 0.1;
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;

points = 4; 
c = linspace(x_b,x_e,points);
Dc = interp1(xsim,Dsim,c,'pchip');

l = ones(points, length(x));
for i = 1:points
    for j = 1:points
        if(i~=j)
            l(i,:) = l(i,:).* (x - c(j))/(c(i)-c(j));
        end 
    end 
end

l = l';

for i = 1:points
    plot(l(:,i));
    hold on;
end