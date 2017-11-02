%%%%%%%%%%%%%%%%%
% load inp_mat; load gamma; load omega_all; dtheta = 1e-8; Bnumber = 1;
% omega = omega_all(1);
%%%%%%%%%%%%%%%%%

function [Jnum] = Test_Jacobian_D(inp_mat,gamma,omega,dtheta,Bnumber)


P = gamma.D;

lri = size(inp_mat.C,1);
Jnum = zeros(lri,length(P));


idxFree = 1:length(gamma.D);

for k=length(idxFree):-1:1
    
    % profielen
    theta = P;
    
    % left side
    thetaL(idxFree(k)) = theta(idxFree(k))*(1-dtheta);
    % right side
    thetaR(idxFree(k))  = theta(idxFree(k))*(1+dtheta);

      
    % thetad (left_side)
    gamma.D = thetaL;
    [Gd,~] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);
    e1 = Gd(:,Bnumber);
    eo1 =[real(e1);imag(e1)];
    

    
    % theta (right_side)
    gamma.D = thetarR;
    [G,~] = Bspline_G_and_J_G_case(inp_mat,gamma,omega);
    e2 = G(:,Bnumber);

    eo2=[real(e2);imag(e2)];
    % Calculate Jacobian
    J(:,idxFree(k))=(eo2-eo1)/(2*dtheta*theta(idxFree(k)));
    Jnum(:,idxFree(k)) = J(1:lri,idxFree(k))+J(lri+1:end,idxFree(k))*1i;
end
% end