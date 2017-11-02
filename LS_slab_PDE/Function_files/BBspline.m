% d = 4; i = 8; t = 0.9; 
% T = [-0.1857, -0.0714, 0.0429, 0.1571, 0.2714, 0.3857, 0.5000,...
% 0.6143, 0.7286, 0.8429, 0.9571, 1.0714, 1.1857];

 function [value] = BBspline(d,i,t,T)
    if(d==0)
        if(t>T(i) && t<=T(i+1))
            value=1;
            return;
        else
            value=0;
            return;
        end;
    end;
    value=(t-T(i))/(T(i+d)-T(i))*BBspline(d-1,i,t,T)+...
        (T(i+d+1)-t)/(T(i+d+1)-T(i+1))*BBspline(d-1,i+1,t,T);
    if(isnan(value))
            value=0;
    end;
    return;
% end
 