function S = supply_function_mb(E,E_FE,E_FC,k_B,theta)
    beta = k_B*theta;
    Col = exp((E_FC-E)/(beta));
    Emi = exp((E_FE-E)/(beta));
    S = ((Emi)-(Col));
end