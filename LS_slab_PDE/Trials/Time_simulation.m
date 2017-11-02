sys = ss(A,B,C,0);
t = linspace(0,400,400);
y = lsim(sys,u,t);
plot(t,y(:,1:end))
legend('y1','y2','y3','y4','y5','y6','y7','y8','y9','y10')