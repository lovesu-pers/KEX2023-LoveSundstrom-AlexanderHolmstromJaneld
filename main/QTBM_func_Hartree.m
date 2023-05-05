% QTBM function for Hartree units. 
% Inputs
% mw: relative effective mass in well region, md: rel.eff.mass in barrier
% d: barrier width, w: well width, xb: start coordinate of barrier
% Nxd: number of points in each segment to calc wave function for
% U0: barrier height, V: bias
% E: incident energy, E0: ref energy
% LE: emitter spacer, LC: collector spacer

% Outputs
% T: transmission coeff. psiLR: out of domain wave function
% xtest: internal x-points, Ux: potential profile
% U_v: internal potential used in T-matrices
% mx_tot: effective mass vector for whole domain, m_v: internal effmass v
% x_vloc: internal vector, 
% psiFWD: wave function moving in positive x-direction
% psiBCK: wave function moving in negative x-direction

function [PSI, T] = QTBM_func_Hartree(xv,U0,V_a,E,xb,w,d,e0,eb,ew,md,mw)
        dx = xv(2)-xv(1);
        Nx = length(xv);
        mL = mw;
        mR = mw;

        m_e = 1;
        a0 = 1; 
        hbar = 1;

        Eh = hbar^2/(m_e*a0^2);


        
        q = Eh;
    
        [Ux, ~, mx, ~, ~, ~, ~,~,~,~] = ... 
                                simpUx_mxV4(xv,U0,V_a,xb,w,d,e0,eb,ew,q,md,mw);

        % Element i system-matris
        a = hbar^2/(4*dx^2);
        s = zeros(Nx+1,1);
        dd = zeros(Nx,1);
        
        for j = 2:Nx-1
            dd(j) = a*( 1/(mx(j-1)) + 2/(mx(j)) + 1/(mx(j+1)) ) + Ux(j);
            s(j) = a*(1/(mx(j-1))+1/(mx(j)));
        end
        
        % Fixartill ränder
        dd(1) = a*( 1/(mL) + 2/(mx(1)) + 1/(mx(2)) ) + Ux(1);
        dd(Nx) = a*( 1/(mx(Nx-1)) + 2/(mx(Nx)) + 1/(mR) ) + Ux(Nx);
        
        s(1) = a*(1/(mL)+1/(mx(1)));
        s(Nx) = a*(1/(mx(Nx))+1/(mx(Nx-1)));
        s(Nx+1) = a*(1/(mR)+1/(mx(Nx)));
        
        % Vågtal och propergator
        k1 = sqrt(2*mx(1)*(E-Ux(1)))/hbar;
        kNx = sqrt(2*mx(Nx)*(E-Ux(Nx)))/hbar;
        
        z1 =  exp(1i*k1*(dx) );
        zNx = exp(1i*kNx*(dx) );
        
        % Element till öppna randvillkor
        alpha1 = 1/(z1 - (z1^(-1)));
        beta1 = (-z1^(-1))/( (z1) - z1^(-1));
        alphaNx = 1/(zNx - (zNx^(-1)));
        betaNx = (zNx^(-1))/( zNx^(-1)-zNx);
        
        % Systemmatris        
        dg1 = [-s(1:Nx);betaNx;0];
        dg2 = [alpha1;dd-E;alphaNx];
        dg3 = [0;beta1;-s(2:Nx+1)];
        Hs3 = spdiags([dg1,dg2,dg3],-1:1,Nx+2,Nx+2);
        
        % Högerled a_1 = 1 inkommande från vänster
        avec = zeros(Nx+2,1);
        avec(1) = 1;
        avec(end) = 0;
        
        % Löserut vågfunktionen
        PSI = Hs3\sparse(avec);
        
        % Amplituder
        an = PSI(Nx);
        a1 = PSI(2)-1;
        
        T = (kNx/k1)*(abs(an))^2;

end