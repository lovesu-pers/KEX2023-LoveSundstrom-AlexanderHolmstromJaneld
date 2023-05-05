function Psi = Rk4_solve(psi0, Ux, hbar, m, k, dx, dt, Nt)
    tic

    Nx = length(psi0);
    
    Psi = zeros(Nx,Nt);
    
    Psi(:,1) = psi0;

    % Kinetisk energi
    T = hbar^2 * k^2 / (2*m);
    
    % Fashastighet vid ränder, behövs till ABC
    vp_1 = (T + Ux(1)) / (sqrt(2*m*T));
    vp_end = (T + Ux(Nx)) / (sqrt(2*m*T));
    
    for j = 1:Nt-1
        psi_t = Psi(:,j);

        k1 = SEdt(psi_t, Ux, hbar, m, dx);
        k2 = SEdt(psi_t+dt*0.5*k1, Ux, hbar, m, dx);
        k3 = SEdt(psi_t+dt*0.5*k2, Ux, hbar, m, dx);
        k4 = SEdt(psi_t+dt*k3, Ux, hbar, m, dx);
        
        K = k1+2*k2+2*k3+k4;

        M1 = psi_t(end); 
        M2 = psi_t(end-1);
        M3 = psi_t(1); 
        M4 = psi_t(2);
        
        Psi(:,j+1) = psi_t + dt*K/6;
        
        ABC_R = (vp_1*dt - dx)/(vp_1*dt + dx);
        ABC_L = (vp_end*dt - dx)/(vp_end*dt + dx);
        %  ABC, absorberande randvillkor
        Psi(Nx,j+1) = M2 + ABC_R * (Psi(Nx-1,j+1) - M1);
        Psi(1,j+1) = M4 + ABC_L * (Psi(2,j+1) - M3);
    end
    toc
end