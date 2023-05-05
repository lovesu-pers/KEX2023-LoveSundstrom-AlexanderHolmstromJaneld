% Crank-Nicolson function for Hartree units. 
% Inputs:
% psi0: inital wave function
% Ux: global potential profile, 
% mx: global mass vector,
% Nx: number of space points
% Nt: total number of time steps
% Nframe: Approximate number of time steps to save if animateflag=1
% dx: grid size in space
% dt: grid size in time
% hbar: set to 1
% animteflag: setting to 1 saves Nframes of the timesteps, 0 only saves
% last frame

% Outputs
% ansmat: Nx by Nframe matrix if animateflag=1, else Nx by 1 vector

function ansmat = CN_solve_mx(psi0, Ux, mx, Nx, Nt, Nframe, dx, dt, hbar,animateflag)
    mL = mx(1);
    mR = mx(end);

    % Element i system-matris
    a = hbar^2/(1*dx^2);
    s = zeros(Nx,1);
    dd = zeros(Nx,1);
        
    for j = 2:Nx-1
        dd(j) = a*(1/(mx(j-1)+mx(j))+1/(mx(j)+mx(j+1)) ) + Ux(j);
        s(j) = a*(1/(mx(j-1)+mx(j)));
    end
    
    % Fixartill ränder
    dd(1) =a*(1/(mL+mx(1))+1/(mx(1)+mx(2)) ) + Ux(1);
    dd(Nx) = a*(1/(mx(Nx-1)+mx(Nx))+1/(mx(Nx)+mR) ) + Ux(Nx);
    
    s(1) = a*(1/(mL+mx(1)));
    s(Nx) = a*(1/(mx(Nx-1)+mx(Nx)));

    A = (spdiags([-[s(2:Nx);0], dd, -s],-1:1,Nx,Nx));
%     A = (spdiags([-s, dd, -s],-1:1,Nx,Nx));
    
    B = speye(Nx,Nx) - (dt/(2i*hbar))*A;
    C = (dt/(2i*hbar))*A + speye(Nx,Nx);
    
    % Svarsmatris och initialvillkor
    ansmat = zeros(Nx,Nframe);
    ansmat(:,1) = psi0;
    psi_p = psi0;
    ff = 2;

    sframe = round(Nt/Nframe);

    % Rumsindelning längs med rader och tidpunkter i kolumner

% tic
    for n = 2:Nt
        psi_p = B\(C*psi_p);
        if mod(n,sframe)==0 && animateflag ==1 
            ansmat(:,ff) = psi_p;
            ff=ff+1
        elseif animateflag == 0 
            ansmat = psi_p;
        end


    end
%     toc
    ansmat(:,end) = psi_p;

end