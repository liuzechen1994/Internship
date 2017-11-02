x_b = 0.1; x_e = 0.9; N = 400;
dx = (x_e-x_b)/(N-1);% N-1 or N error has serious consequences needs to be investigated
x = x_b:dx:x_e;
degree = 7;

T  = zeros(length(x),degree);
% scale -1 to 1 
x_scale= 2/(x(end) - x(1)) * ( x - 0.5);

T(:,1) = 1;
T(:,2) = x_scale;

for i = 3:degree
    T(:,i) = 2 * x_scale' .*T(:,i-1) - T(:,i-2);

end

plot(x,T)
legend('1','2','3','4','5','6','7')

% fplot(chebyshevT(0:6, x))
% axis([0.1 0.9 -2 2])
% grid on
% 
% ylabel('T_n(x)')
% legend('T_0(x)', 'T_1(x)', 'T_2(x)', 'T_3(x)', 'T_4(x)',  'T_5(x)', 'T_6(x)','Location', 'Best')
% title('Chebyshev polynomials of the first kind')