function T = T_prim_linear(x, beta, alpha, m, E, hbar)
u = ( (2*alpha*m/(hbar.^2)).^(1/3)*( alpha + beta*x - E) )/beta;
alpha + beta*x - E
du_i_dx = (2*alpha*m/(hbar.^2)).^(1/3);

psi_lin_A = airy(0,u);
dpsi_lin_A_dx = du_i_dx*airy(1,u);
psi_lin_B = airy(2,u);
dpsi_lin_B_dx = du_i_dx*airy(3,u);

T = [ psi_lin_A, psi_lin_B; dpsi_lin_A_dx, dpsi_lin_B_dx ];

end

