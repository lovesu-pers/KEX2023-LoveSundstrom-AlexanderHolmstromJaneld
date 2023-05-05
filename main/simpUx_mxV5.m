% Creates RTD potential profile, effective mass vector, and relevant index 
% Inputs: 
% xv: space vector, 
% U_0 = conduction band offset / barrier height
% x_b: coordinate for start of first barrier,
% w: well width, d: barrier width
% LE: Emitter spacer length, LC: Collector spacer length
% q: fundamental charge, set to 1 for Hartree units,
% mw: relative effective mass in well region, md: rel.eff.mass in barrier
% E_0: Base energy, can be set to 0 or bottom of C-band

% Output: 
% Ux: global potential profile, Ux_loc: local potential profile
% mx: global mass vector, mx_loc: local mass vector,
% NL: number of points left of local system, 
% NR: number of points right of local system

function [Ux, Ux_loc, mx, mx_loc, x_loc, NL, NR,IDd1,IDw,IDd2] = ... 
                            simpUx_mxV5(xv,U_0,V_a,x_b,w,d,LE,LC,q,md,mw,E_0)
    Nx = length(xv);
    x0 = xv(1);
    xfin = xv(end);

    L = xfin-x0;
    L2 = 2*d+w+LE+LC;
    dx = xv(2)-xv(1);
    
    Ux = zeros(Nx, 1);
    
    [~,LEstart] = min(abs(xv-(x_b-LE)));
    [~,d1start] = min(abs(xv-(x_b)));
    LEfin = d1start-1;

    [~,d1fin] = min(abs(xv-(x_b+d)));
    wstart = d1fin+1;
    [~,wfin] = min(abs(xv-(x_b+w+d)));
    wfin = wfin-1;
    d2start = wfin+1;
    [~,d2fin] = min(abs(xv-(x_b+w+2*d)));
    LCstart = d2fin+1;
    [~,LCfin] = min(abs(xv-(x_b+w+2*d+LC)));
    LCfin = LCfin-1;
    
    IDd1 = [d1start, d1fin];
    IDw = [wstart,wfin];
    IDd2 = [d2start,d2fin];


    xb1 = round([(x_b-x0)/dx + 1, (x_b+d-x0)/dx + 1]);
    xw  = round([xb1(2)+1, xb1(2)+w/dx+1]);
    xb2 = round([xw(2)+1, xw(2) + d/dx+1]);
    

    Ux_b1 = d1start:d1fin;
    Ux_w = wstart:wfin; 
    Ux_b2 = d2start:d2fin;
    
    
    Ux(1:LEstart-1) = 0;
    Ux(Ux_b1) = U_0 ;
    Ux(Ux_b2) = U_0;
    Ux(Ux_w) =   0;
    Ux(LCfin+1:end) = -V_a*q; % Varning
    
    UV = zeros(Nx,1);
    UV(LEstart:LCfin) = linspace(0,-V_a*q,LCfin-LEstart+1).';
    Ux = Ux+UV+E_0;

    Ux_loc = Ux(LEstart-1:LCfin+1);
    x_loc = xv(LEstart-1:LCfin+1);


    [~,d1start_loc] = min(abs(x_loc-(x_b)));
    [~,d1fin_loc] = min(abs(x_loc-(x_b+d)));
    wstart_loc = d1fin_loc+1;
    [~,wfin_loc] = min(abs(x_loc-(x_b+w+d)));
    d2start_loc = wfin_loc+1;
    [~,d2fin_loc] = min(abs(x_loc-(x_b+w+2*d)));

    N = length(Ux_loc)-2;
    mx_loc = mw*ones(N+2,1);
    mx_loc(1) = mw;
    mx_loc(end) = mw;
    mx_loc(d1start_loc-1:d1fin_loc-1) = md;
    mx_loc(d2start_loc-2:d2fin_loc-1) = md;
    
    NL = LEstart-1;
    NR = Nx-(NL+length(Ux_loc)); 
    % NR2 = IDd
    mxL = mw*ones(NL, 1);
    mxR = mw*ones(NR, 1);

    mx = [mxL;mx_loc;mxR];
end