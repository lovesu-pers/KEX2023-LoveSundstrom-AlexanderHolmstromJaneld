% QTBM function for Hartree units. 
% Inputs
% xv: space vector
% U_0: conduction band offset / barrier height
% Vb: bias voltage
% E: incident energy
% xb: start coordinate of barrier
% d: barrier width, w: well width,
% LE: emitter spacer, LC: collector spacer
% mw: relative effective mass in well region, md: rel.eff.mass in barrier
% Nxd: number of points in each segment to calc wave function for
% E: incident energy, 
% E0:  ref energy, can be set to 0 or bottom of C-band
% aE: set to 1 for emitter incident (from right), 0 for left/conduction
% aC: set to 1 for collector incident (from left), 0 for right/emitter
% Vflag: set to zero,
% Vtot: can be left empty


% Outputs
% PSI: wave function
% Transmission coefficient

function [PSI, T] = QTBM_func_Hartree_V2(xv,U0,Vb,E,xb,w,d,LE,LC,...
    md,mw,E0,aE,aC,Vflag,Vtot)

        dx = xv(2)-xv(1);
        Nx = length(xv)-2;
        mL = mw;
        mR = mw;

        m_e = 1;
        a0 = 1; 
        hbar = 1;

        Eh = hbar^2/(m_e*a0^2);


        
        q = Eh;
        
        if Vflag == 1
            Ux = Vtot;
            [~, ~, mx, ~, ~, ~, ~,~,~,~] = ... 
                                simpUx_mxV5(xv,U0,Vb,xb,w,d,LE,LC,q,md,mw,E0);
        else
        [Ux, ~, mx, ~, ~, ~, ~,~,~,~] = ... 
                                simpUx_mxV5(xv,U0,Vb,xb,w,d,LE,LC,q,md,mw,E0);
%         Ux = Ux+E_0;
        end

%         if aC==1 && aE==0
%             [Ux, ~, mx, ~, ~, ~, ~,~,~,~] = ... 
%                                 simpUx_mxV5(xv,U0,V_a,xb,w,d,LE,LC,q,md,mw,E_b);
%             Ux = Ux+V_a*q;
%             Ux = flip(Ux);
%             mx = flip(mx);
%             aC = 0;
%             aE = 1;
%              
%              E = (E+V_a*q);
%         end
    
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
        if aC == 1 && E<Ux(1)
%             kNx=-kNx;
        end
        
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
        avec(1) = aE;
        avec(end) = -aC;
        
        % Löserut vågfunktionen
        PSI = Hs3\sparse(avec);
        
        % Amplituder
        if aE==1
            an = PSI(Nx);
            
            T = (abs(kNx)/abs(k1))*(abs(an))^2;
        else
            an = PSI(1);
            a1 = PSI(Nx);
            T = (abs(k1)/abs(kNx))*(abs(an))^2;
        end
        

end