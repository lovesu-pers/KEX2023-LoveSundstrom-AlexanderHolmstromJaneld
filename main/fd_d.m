function f = fd_d(E,Ef,kB,T)
    f = 1./(1+exp((E-Ef)/(kB*T)));
end