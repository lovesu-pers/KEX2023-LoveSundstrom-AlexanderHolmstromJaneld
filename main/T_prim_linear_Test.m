function T = T_prim_linear_Test(x, Vcurr, Vnext, xcurr, xnext, m, E, hbar)
%% Wolfram Alpha, with V = kx + m
% alpha = Vcurr-xcurr*((Vnext-Vcurr)/(xnext-xcurr));
% beta = (Vnext-Vcurr)/(xnext-xcurr);
% 
% u = ( (2*alpha*m/(hbar.^2)).^(1/3)*( alpha + beta*x - E) )/beta;
% %alpha + beta*x - E
% du_i_dx = (2*alpha*m/(hbar.^2)).^(1/3);
% 
% psi_lin_A = airy(0,u);
% dpsi_lin_A_dx = du_i_dx*airy(1,u);
% psi_lin_B = airy(2,u);
% dpsi_lin_B_dx = du_i_dx*airy(3,u);
% 
% T = [ psi_lin_A, psi_lin_B; dpsi_lin_A_dx, dpsi_lin_B_dx ];
%% Wayne W. Lui and Masao Fukuma
% F = -(Vnext - Vcurr)/(xnext-xcurr);
% r = -(2*m*F/(hbar.^2))^(1/3);
% c = xcurr + (Vcurr-E)/F;
% z = r*(c-x);
% zprim = -r;
%% Gehring Andreas
% V = Vcurr + (x-xcurr)*( (Vnext-Vcurr)/(xnext-xcurr) );
% z = -(2*m/((hbar.^2))).^(1/3)*( (xnext-xcurr)/(Vnext-Vcurr) )^(2/3)*( E - V );
% zprim = -(2*m/((hbar.^2)))^(1/3)*( (xnext-xcurr)/(Vnext-Vcurr) )^(1/3);
%% Wolfram Alpha
% z = (2).^(1/3)*(m*(Vcurr-Vnext)/(hbar*(xcurr-xnext))).^(1/3)*(Vnext*(xcurr-x) + Vcurr*(x-xnext) + E*(xnext-xcurr))/(Vcurr-Vnext);
% zprim = (2).^(1/3)*(m*(Vcurr-Vnext)/(hbar*(xcurr-xnext)));

%% Christian Jirauschek, Accuracy of Transfer Matrix Approaches for Solving the Effective Mass Schrödinger Equation
Vzj = (Vnext-Vcurr)/(xnext-xcurr);
epsilon = ( hbar.^2*Vzj.^2/(2*m) ).^(1/3);
s = (Vcurr-E)/epsilon;
l = epsilon/Vzj;

z = s + (x-xcurr)/l;
zprim = 1/l;

T = [ airy(0,z), airy(2,z); 1/m*zprim*airy(1,z), 1/m*zprim*airy(3,z) ];


end

