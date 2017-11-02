g = sdpvar(7,1);
L = [inp_mat.BP.P * g >= 0];
optimize(L);
g = value(g)