function ansmat = CN_solve(psi0, Ux, Nx, Nt, dx, dt, m)


    hbar = 1;
    a1 = hbar/(m*dx^2) * ones(Nx+1,1);
    a2 = -hbar/(2*m*dx^2) * ones(Nx+1,1);
    a3 = spdiags(Ux,0,Nx+1,Nx+1);

    A = spdiags([a2, a1, a2],-1:1,Nx+1,Nx+1)  ...
        + a3;
    
    % Om man vill ha periodiska randvillkor
    A(end,1) = -hbar/(2*m*dx^2);
    A(1,end) = -hbar/(2*m*dx^2);
    
    B = speye(Nx+1,Nx+1) - (dt/(2i*hbar))*A;
    C = (dt/(2i*hbar))*A + speye(Nx+1,Nx+1);
    
    % Svarsmatris och initialvillkor
    ansmat = zeros(Nx+1,Nt);
    ansmat(:,1) = psi0;
    
    % Rumsindelning längs med rader och tidpunkter i kolumner
    tic
    for n = 2:Nt
        ansmat(:,n) = B\(C*ansmat(:,n-1));
    end
    toc

end