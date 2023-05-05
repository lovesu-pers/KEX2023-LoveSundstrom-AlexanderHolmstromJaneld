clear
close all
clc
% Program för att beräkna transmission koefficient för en dubbelbarriär
% genom transfer matris metod vid antagande av en lokalt linjär potential
% inuti barriären och stående vågor utanför
%%

% eV = 1;
% m_e = 1;
% nm = 1;
% fs = 1;
% hbar = 0.658217*eV*fs;
% 
% mw = 1; %0.067*m_e;
% md = 1; %(0.063+0.083*0.3)*m_e; %0.067*m_e;
% 
% 
% U0 = 0.38*eV  % Refhöjd
% V = 0.5;
% E = 0.25
% 
% d = 2.5*nm;     % Barriärtjocklek
% w = 5*nm;       % Brunnbredd
% L = 2*d+w;      % Längd på barriär+brunn område
% x_b = 0;        % Start x-koordinat av barriär 1
% 
% Nbias = 50;
% Biasv = linspace(0.1,1,Nbias);
% 
% Nx = 1000;  % Antal interna segment
% xv = linspace(-3*L,3*L,Nx)';
% dx = xv(2)-xv(1);



m_e = 1;
a0 = 1; 
hbar = 1;

Eh = hbar^2/(m_e*a0^2);

fs = hbar/Eh; 
kB =  3.167e-6; % Eh/K
q = Eh;
A0 = q*m_e*kB^2/(2*pi^2*hbar^3); % Richardson konstant A/(m^2 K^2)
e0 = 1;
eb= 1;
ew = 1;

% Hartree enheter energi=Eh, längd=a_0=Bohr-radie osv
q_SI = 1.602176634e-19;       % C
a0_SI = 5.291772109e-11;      % m
Eh_SI = 4.35974472220712e-18; % J
hbar_SI = 1.054571817e-34;    % Js
m_e_SI = 9.109534e-31;        % kg
V_Hartree_SI = 27.2113862456; % V
nm_a0_sf = 1e-9/a0_SI ;
j_sf = q_SI*Eh_SI/(hbar_SI*a0_SI^2);

T = 300; % Kelvin
EF = 0;  % Ferminivå  

factor = 1e-1;

% Konstanter
mw = 0.067*m_e;
md = 0.0919*m_e;        %(0.063+0.083*0.3)*m_e; %0.067*m_e;
d = 3*nm_a0_sf*a0*factor;     % Barriärtjocklek
w = 4*nm_a0_sf*a0*factor; %4*nm_a0_sf*a0;       % Brunnbredd
L = (2*d+w);      % Längd på barriär+brunn område
xb = 0;
mL = 1*mw;
mR = 1*mw;
Amod =  (mw*A0*T/(kB*m_e) );

% Antal spänningspunkter
Nbias = 50;
NE = 200;

% Barriärhöjd = 0.2 eV = 0.007 Eh
U0 = (0.231*q_SI/Eh_SI)*Eh;

% Partikelenergi
Ev = linspace(0.01,0.5,NE)*q_SI/Eh_SI;
E = Ev(1)

%Biasv = linspace(1/Nbias,1+1/Nbias,Nbias)/V_Hartree_SI;
Biasv = linspace(0,1,Nbias)/V_Hartree_SI;

Nx = 1000;  % Antal interna segment
xv = linspace(-3*L,3*L,Nx)';
dx = xv(2)-xv(1);

% %Potential (global Ux och lokal U_v, mx m_v pss)
% [Ux,U_v,mx_tot,m_v, x_vloc,NL,NR, sb1, sb2, sw,d1start,wstart,d2start,d2fin] =simpUx_mx(xv,U0,V,0,w,d,1,1,1,1,md,mw);
% figure
% plot(U_v)
% figure
% plot(xv,Ux)

% Kollar potentialen med max framspänning
% [Ux, Ux_loc, mx, mx_loc, x_loc, NL, NR,IDd1,IDw,IDd2] = ... 
%                         simpUx_mxV4(xv,U0,Biasv(end),xb,w,d,1,1,1,q,md,mw);
% 
% figure(1)
% plot(xv, Ux);
% yline(E)

% Upplösning i x osv
dx = 0.1;
x0 = -5; 
xfin = L+5; 
Nx = round((xfin-x0)/dx);   % Antal x-punkter till QTBM
dx = (xfin-x0) / Nx;
xv = x0:dx:xfin-dx;
NxTMM = 2000;



%%

Transvec_analytical = zeros(Nbias,1);
Transvec_TMM = zeros(Nbias,1);

TvM = zeros(Nbias,1);

for element=20:20
    V = Biasv(element);
    format 
   % Potential (global Ux och lokal U_v, mx m_v pss)
    [Ux,U_v,mx_tot,m_v, x_vloc,NL,NR, sb1, sb2, sw,d1start,wstart,d2start,d2fin] =simpUx_mx(xv,U0,V,0,w,d,1,1,1,1,md,mw);
    plot(xv,Ux)
    % Övergång mellan emitter och dubbelbarriär
    xstartb1 = xv(d1start)
    xendb1 = xv(wstart-1)
    V_at_xstartb1 = Ux(d1start)
    V_at_xendb1 = Ux(wstart-1)
    %beta3 = -sb1;
    %alpha3 = V_at_xstartb1-beta3*xstartb1;
    emitterV = U_v(1)
    kL = sqrt(2*mw*( E-emitterV ))/hbar;

    if V_at_xstartb1 == V_at_xendb1
        T_prim_const_3_atd1start = T_prim_const(xstartb1,V_at_xstartb1,md,E,hbar,1);
        T_prim_const_L_atd1start = T_prim_const(xstartb1,emitterV,mw,E,hbar,1);
        T_L = T_prim_const_L_atd1start\T_prim_const_3_atd1start;
        %T_L = T_prim_const_3_atd1start\T_prim_const_L_atd1start;
    else
    
        T_prim_lin_3_atd1start = T_prim_linear_Test(xstartb1, V_at_xstartb1, V_at_xendb1, xstartb1, xendb1, md, E, hbar);
        T_prim_const_L_atd1start = T_prim_const(xstartb1,emitterV,mw,E,hbar,0);
        T_L = T_prim_const_L_atd1start\T_prim_lin_3_atd1start;
        %T_L = T_prim_lin_3_atd1start\T_prim_const_L_atd1start;
    end
    % Övergång mellan första barriär och brunn
    xstartw = xv(wstart)
    xendw = xv(d2start-1)
    V_at_xstartw = Ux(wstart)
    V_at_xendw = Ux(d2start-1) 
    %beta2 = sw;
    %alpha2 = V_at_xstartw-beta2*xstartw;
    
    if V_at_xstartb1 == V_at_xendb1
        T_prim_const_2_atd1end = T_prim_const(xendb1,V_at_xstartw,mw,E,hbar,1);
        T_prim_const_3_atd1end = T_prim_const(xendb1,V_at_xstartb1,md,E,hbar,1);
        
        T_2 = T_prim_const_3_atd1end\T_prim_const_2_atd1end;
        %T_2 = T_prim_const_2_atd1end\T_prim_const_3_atd1end;
    else
    
        T_prim_lin_3_atd1end = T_prim_linear_Test(xendb1, V_at_xstartb1, V_at_xendb1, xstartb1, xendb1, md, E, hbar);
        T_prim_lin_2_atd1end = T_prim_linear_Test(xendb1, V_at_xstartw, V_at_xendw, xstartw, xendw, mw, E, hbar);
        
        T_2 = T_prim_lin_3_atd1end\T_prim_lin_2_atd1end;
        %T_2 = T_prim_lin_2_atd1end\T_prim_lin_3_atd1end;
    end
    
    % Övergång mellan brunn och andra barriär
    xstartd2 = xv(d2start)
    xendd2 = xv(d2fin)
    V_at_xstartd2 = Ux(d2start)
    V_at_xendd2 = Ux(d2fin)
    %beta = -sb2;
    %alpha1 = V_at_xstartb2-beta*xstartb2;
    if V_at_xstartb1 == V_at_xendb1
        T_prim_const_1_atd2start = T_prim_const(xstartd2,V_at_xstartd2,md,E,hbar,1);
        T_prim_const_2_atd2start = T_prim_const(xstartd2,V_at_xstartw,mw,E,hbar,1);
        
        T_1 = T_prim_const_2_atd2start\T_prim_const_1_atd2start;
        %T_1 = T_prim_const_1_atd2start\T_prim_const_2_atd2start;
    else
        T_prim_lin_2_atd2start = T_prim_linear_Test(xstartd2, V_at_xstartw, V_at_xendw, xstartw, xendw, mw, E, hbar);
        T_prim_lin_1_atd2start = T_prim_linear_Test(xstartd2, V_at_xstartd2, V_at_xendd2, xstartd2, xendd2, md, E, hbar);
        
        T_1 = T_prim_lin_2_atd2start\T_prim_lin_1_atd2start;
        %T_1 = T_prim_lin_1_atd2start\T_prim_lin_2_atd2start;
    end

    % Övergång mellan dubbelbarriär och collector

    collectorV = U_v(end)
    kR = sqrt(2*mw*( E-collectorV ))/hbar;
    if V_at_xstartb1 == V_at_xendb1
        T_prim_const_1_atd2end = T_prim_const(xendd2,V_at_xstartd2,md,E,hbar,1);
        T_prim_const_R_atd2end = T_prim_const(xendd2,collectorV,mw,E,hbar,1);
        T_R = T_prim_const_1_atd2end\T_prim_const_R_atd2end;
        %T_R = T_prim_const_R_atd2end\T_prim_const_1_atd2end;
         
    else
        T_prim_lin_1_atd2end = T_prim_linear_Test(xendd2, V_at_xstartd2, V_at_xendd2, xstartd2, xendd2, md, E, hbar);
        T_prim_const_R_atd2end = T_prim_const(xendd2,collectorV,mw,E,hbar,0);
        T_R = T_prim_lin_1_atd2end\T_prim_const_R_atd2end;
        %T_R = T_prim_const_R_atd2end\T_prim_lin_1_atd2end;
    end
    % Totala transfer matris
    T = T_L*T_2*T_1*T_R;
    %T = T_R*T_1*T_2*T_L;
    
    Transfer_coefficient = abs(kR)/(abs(kL))*1/(abs(T(1,1))^2);
    %Transfer_coefficient = abs(kR)/(abs(kL))*(abs(T(1,1)-T(1,2)*T(2,1)/T(2,2))^2);

    Transvec_TMM(element) = Transfer_coefficient;

%     % test enkel barriär
%     p = 1e-0;
%     Vstart = 0;
%     V1 = 2*p;
%     V2 = 1.5*p;
%     Vend = 0;
%     xstart = 0;
%     xend = 1*p;
%     E = Espec(element)*p;
%     % Övergång 1
%     
%     T_const_L_1 = T_prim_const(xstart,Vstart,m_e,E,hbar);
%     T_const_B_1 = T_prim_const(xstart,V1,m_e,E,hbar);
%     %T_lin_L = T_prim_linear_Test(xstart, V1, V2, xstart, xend, m_e, E, hbar);
%     T_L = inv(T_const_L_1)*T_const_B_1;
%     
%     % Övergång 2
%     
%     T_const_R_2 = T_prim_const(xend,Vend,m_e,E,hbar);
%     T_const_B_2 = T_prim_const(xend,V1,m_e,E,hbar);
%     %T_lin_R = T_prim_linear_Test(xend, V1, V2, xstart, xend, m_e, E, hbar);
%     T_R = inv(T_const_B_2)*T_const_R_2;
% 
%     T = T_L*T_R;
%     
%     Transfer_coefficient_TMM = 1/(abs(T(1,1))^2);
%     Transvec_TMM(element) = Transfer_coefficient_TMM;
% 
%     % Analytical
%     a = xend-xstart;
%     k0 = sqrt(2*m_e*E)./hbar;
%     k1 = sqrt(2*m_e*(E-V1))./hbar;
%     Transfer_coefficient_analytical = abs(4*k0*k1*exp(-1i*a*(k0-k1)) /( (k0+k1).^2 - exp(2i*a*k1)*(k0-k1).^2 )).^2;
%     Transvec_analytical(element) = Transfer_coefficient_analytical;
% 



%     V = Biasv(element);
%      
%     % Potential (global Ux och lokal U_v, mx m_v pss)
%     [Ux,U_v,mx_tot,m_v, x_vloc,NL,NR, sb1, sb2, sw,d1start,wstart,d2start,d2fin] =simpUx_mx(xv,U0,V,0,w,d,1,1,1,1,md,mw);
%     plot(xv,Ux)
%     %plot(xv,Ux)
% 
%     Vemitter = U_v(1);
%     Vcollector = U_v(end);
% 
%     m_L = m_e;
%     m_R = m_e;
%     
%     kL = sqrt(2*m_L*(E-Vemitter))/hbar;
%     kR = sqrt(2*m_R*(E-Vcollector))/hbar;
%     
%     xstartd1 = xv(d1start);
%     xendd1 = xv(wstart-1);
%     Vstartd1 = Ux(d1start);
%     Vendd1 = Ux(wstart-1);
%     m_3 = m_e;
% 
%     xstartw = xv(wstart);
%     xendw = xv(d2start-1);
%     Vstartw = Ux(wstart);
%     Vendw = Ux(d2start-1);
%     m_2 = m_e;
%    
%     xstartd2 = xv(d2start);
%     xendd2 = xv(d2fin);
%     Vstartd2 = Ux(d2start);
%     Vendd2 = Ux(d2fin);
%     m_1 = m_e;
% 
%     x3 = xv(d1start);
%     x2 = xv(wstart);
%     x1 = xv(d2start);
%     x0 = xv(d2fin);
%   
%     % övergång mellan emittor och första barriär, i x = x3
%     % Bidrag från tillstånd 3
%     V3_x3 = Vstartd1 + (x3-xstartd1)*( (Vendd1-Vstartd1)/(xendd1-xstartd1) );
%     z3_x3 = -(2*m_3/((hbar.^2)))*( (xendd1-xstartd1)/(Vendd1-Vstartd1) )^(2/3)*( E - V3_x3 );
%     zprim3_x3 = -(2*m_3/((hbar.^2)))*( (xendd1-xstartd1)/(Vendd1-Vstartd1) )^(1/3);
%     
%     % övergång mellan första barriär och brunn, i x = x2
%     % Bidrag från tillstånd 3
%     V3_x2 = Vstartd1 + (x2-xstartd1)*( (Vendd1-Vstartd1)/(xendd1-xstartd1) );
%     z3_x2 = -(2*m_3/((hbar.^2)))*( (xendd1-xstartd1)/(Vendd1-Vstartd1) )^(2/3)*( E - V3_x2 );
%     zprim3_x2 = -(2*m_3/((hbar.^2)))*( (xendd1-xstartd1)/(Vendd1-Vstartd1) )^(1/3);
%     
%     % Bidrag från tillstånd 2
%     V2_x2 = Vstartw + (x2-xstartw)*( (Vendw-Vstartw)/(xendw-xstartw) );
%     z2_x2 = -(2*m_3/((hbar.^2)))*( (xendw-xstartw)/(Vendw-Vstartw) )^(2/3)*( E - V3_x2 );
%     zprim2_x2 = -(2*m_3/((hbar.^2)))*( (xendw-xstartw)/(Vendw-Vstartw) )^(1/3);
% 
% 
%     % övergång mellan första brunn och andra barriär, i x = x1
%     % Bidrag från tillstånd 2
%     V2_x1 = Vstartw + (x1-xstartw)*( (Vendw-Vstartw)/(xendw-xstartw) );
%     z2_x1 = -(2*m_3/((hbar.^2)))*( (xendw-xstartw)/(Vendw-Vstartw) )^(2/3)*( E - V3_x2 );
%     zprim2_x1 = -(2*m_3/((hbar.^2)))*( (xendw-xstartw)/(Vendw-Vstartw) )^(1/3);
% 
%     
%     % Bidrag från tillstånd 1
%     V1_x1 = Vstartd2 + (x1-xstartd2)*( (Vendd2-Vstartd2)/(xendd2-xstartd2) );
%     z1_x1 = -(2*m_3/((hbar.^2)))*( (xendd2-xstartd2)/(Vendd2-Vstartd2) )^(2/3)*( E - V3_x2 );
%     zprim1_x1 = -(2*m_3/((hbar.^2)))*( (xendd2-xstartd2)/(Vendd2-Vstartd2) )^(1/3);
% 
%     % övergång mellan andra barriär och collector, i x = x0
%     % Bidrag från tillstånd 1
%     V1_x0 = Vstartd2 + (x0-xstartd2)*( (Vendd2-Vstartd2)/(xendd2-xstartd2) );
%     z1_x0 = -(2*m_3/((hbar.^2)))*( (xendd2-xstartd2)/(Vendd2-Vstartd2) )^(2/3)*( E - V3_x2 );
%     zprim1_x0 = -(2*m_3/((hbar.^2)))*( (xendd2-xstartd2)/(Vendd2-Vstartd2) )^(1/3);
% 
%     TL_11 = 1/(2i*kL)*( 1i*kL*airy(0,z3_x3) + zprim3_x3*airy(1,z3_x3) )*exp(-1i*kL*x3);
%     T2_11 = pi*( airy(3,z3_x2)*airy(0,z2_x2) - zprim2_x2/zprim3_x2*airy(1,z2_x2)*airy(2,z3_x2) );
%     TL_12 = 1/(2i*kL)*( 1i*kL*airy(2,z3_x3) + zprim3_x3*airy(3,z3_x3) )*exp(-1i*kL*x3);
%     T2_21 = pi*( -airy(1,z3_x2)*airy(z2_x2) + zprim2_x2/zprim3_x2*airy(1,z2_x2)*airy(0,z3_x2) );
% 
%     T1_11 = pi*( airy(3,z2_x1)*airy(0,z1_x1) - zprim1_x1/zprim2_x1*airy(1,z1_x1)*airy(2,z2_x1) );
% 
%     T2_12 = pi*( airy(3,z3_x2)*airy(2,z2_x2) - zprim2_x2/zprim3_x2*airy(3,z2_x2)*airy(2,z3_x2) );
% 
%     T2_22 = pi*( -airy(1,z3_x2)*airy(2,z2_x2) + zprim2_x2/zprim3_x2*airy(3,z2_x2)*airy(0,z3_x2));
% 
%     T1_21 = pi*( -airy(1,z2_x1)*airy(0,z1_x1) + zprim1_x1/zprim2_x1*airy(1,z1_x1)*airy(0,z2_x1) );
% 
%     TR_11 = pi/(zprim1_x0)*( zprim1_x0*airy(3,z1_x0) - 1i*kR*airy(2,z1_x0) )*exp(1i*kR*x0);
% 
%     T1_12 = pi*( airy(3,z2_x1)*airy(2,z1_x1) - zprim1_x1/zprim2_x1*airy(3,z1_x1)*airy(2,z2_x1) );
% 
%     T1_22 = pi*( -airy(1,z2_x1)*airy(2,z1_x1) + zprim1_x1/zprim2_x1*airy(3,z1_x1)*airy(0,z2_x1) );
% 
%     TR_21 = pi/(zprim1_x0)*( -zprim1_x0*airy(1,z1_x0) + 1i*kR*airy(0,z1_x0) );
% 
%      T11 = ( (TL_11*T2_11 + TL_12*T2_21)*T1_11 + (TL_11*T2_12 + TL_12*T2_22)*T1_21)*TR_11 ...
%         + ( (TL_11*T2_11 + TL_12*T2_21)*T1_12 + (TL_11*T2_12 + TL_12*T2_22)*T1_22)*TR_21;
% 
%      Transfer_coefficient = abs(kR)/(abs(kL))*1/(abs(T11)^2);
%      Transvec_TMM(element) = Transfer_coefficient;
% 
%     zvec = [abs(z3_x3), abs(z3_x2), abs(z2_x2), abs(z2_x1), abs(z1_x1)];
%     zmax = max(zvec);
%     
%     element
%     airy(0,zmax)
%     airy(1,zmax)
%     airy(2,zmax)
%     airy(3,zmax)
end


% figure
% semilogy(Biasv*V_Hartree_SI, Transvec_TMM)
% hold on
%semilogy(Biasv*V_Hartree_SI, TvM)
% hold on
% plot(Vspec,Transvec_analytical)
% legend('TMM', 'Analytical')