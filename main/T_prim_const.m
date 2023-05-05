function T = T_prim_const(x, V, m, E, hbar,status)
k = sqrt(2*m*( E-V ))/hbar;
if status == 1
    T = [exp(1i*k*x), exp(-1i*k*x); ...
        1/m*1i*k*exp(1i*k*x), 1/m*(-1i*k*exp(-1i*k*x))
        ];
else
    T = [exp(1i*k*x), exp(-1i*k*x); ...
        1/m*1i*k*exp(1i*k*x), 1/m*(-1i*k*exp(-1i*k*x))
        ];
end
end

