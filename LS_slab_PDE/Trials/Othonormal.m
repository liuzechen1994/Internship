syms x
fun = @(x) (sin(6*x)/norm(sin(6*x))).*(cos(1*x)/norm(cos(3*x)));
integral(fun,0,1)