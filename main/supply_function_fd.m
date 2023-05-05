function S = supply_function_fd(E,E_FE,E_FC,k_B,theta)
    beta = k_B*theta;
    Col = exp((E_FC-E)/(beta));
    Emi = exp((E_FE-E)/(beta));
    S = log((1+Emi)./(1+Col));
end