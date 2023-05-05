% TMMC function for Hartree units. 
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

function [T,psiLR,xtest,Ux,U_v,mx_tot,m_v, x_vloc,psiFWD,psiBCK] = ...
    TMM_func_Hartree_V2(mw,md,d,w,xb,xv,Nxd,U0,V,E,E0,LE,LC) 

    % Settings for Hartree atomic units
    hbar = 1;
    q = 1;
  
    
    
    [Ux, U_v, mx_tot, m_v, x_vloc, NL, NR,IDd1,IDw,IDd2] = ... 
                            simpUx_mxV5(xv,U0,V,xb,w,d,LE,LC,q,md,mw,E0);

    N = length(U_v)-2;
    xL = xv(1:NL);
    xR = xv(end-NR:end);
    dx = xv(2)-xv(1);
                  
    
    
    psiLR = zeros(length(xL)+length(xR)+N*Nxd);
     
    
    k_v = zeros(N+2,1);
    Delta_v = zeros(N+2,1);
    D_v = zeros(2,2,N+2);
    P_v = zeros(2,2,N+2);
    
    d_v = dx*ones(N+2,1);
    d_v(1) = 0;d_v(end) = 0;
    
    
    
    
    xI = zeros(N+1,1);
    xI(end) = sum(d_v)-d_v(2);
    for xm = 2:N+1
        xI(xm) = sum(d_v(2:xm));
    end
    
    for kj=1:N+2
        k_v(kj) = sqrt(2*m_v(kj)*(E-U_v(kj)))/hbar;
    end
    
    
    for dj = 1:N+1
        Delta_v(dj) = (m_v(dj)/m_v(dj+1))*(k_v(dj+1)/k_v(dj));
        D_v(:,:,dj) = 0.5*[1+Delta_v(dj), 1-Delta_v(dj);
                           1-Delta_v(dj), 1+Delta_v(dj)];
    end
    
    for pj = 2:N+1
        P_v(:,:,pj) = [exp(-1i*k_v(pj)*d_v(pj)), 0; 
                   0,  exp(1i*k_v(pj)*d_v(pj))];  
    end
    
    
    
    T_tot = D_v(:,:,1);
    for Tj = 2:N+1
        T_tot = T_tot*P_v(:,:,Tj)*D_v(:,:,Tj);
    end
    A = 1;
    F = A/T_tot(1,1);
    B = T_tot(2,1)*F;
    
    psiL = A*exp(1i*k_v(1)*xL) + B*exp(-1i*k_v(1)*xL);
    psiR = F*exp(1i*k_v(end)*xR);

    psifwdL = A*exp(1i*k_v(1)*xL) ;
    psibckL = B*exp(-1i*k_v(1)*xL) ;
    psifwdR = F*exp(1i*k_v(end)*xR) ;
    psibckR = zeros(length(psiR),1);
    
    
    
    C = [F;0];
    psitemp = [];
    psifwd = [];
    psibck = [];
    for mm = flip(2:N+1)
        xtemp = linspace(xI(mm-1), xI(mm), Nxd+1).';
        xtemp = xtemp(2:end);
        C = P_v(:,:,mm)*D_v(:,:,mm)*C;
        psiii = C(1)*exp(1i*k_v(mm)*(xtemp-xI(mm-1)))+ ... 
            C(2)*exp(-1i*k_v(mm)*(xtemp-xI(mm-1)));
        psifwd = [C(1)*exp(1i*k_v(mm)*(xtemp-xI(mm-1)));psifwd];
        psibck = [C(2)*exp(-1i*k_v(mm)*(xtemp-xI(mm-1)));psibck];
        psitemp = [psiii; psitemp];
    end
    
    
    
    psiLR = [psiL;psitemp;psiR]; 
    psiFWD = [psifwdL;psifwd;psifwdR];
    psiBCK = [psibckL;psibck;psibckR];

    
    T = 1-( abs(T_tot(2,1))^2 / abs(T_tot(1,1)).^2 );
    
    
    
    x_intrn = linspace(xI(1), xI(end), N*Nxd)';
    xtest = [xL;x_intrn;xR];

end