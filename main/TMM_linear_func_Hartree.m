function Transfer_coefficient = TMM_linear_func_Hartree(mw,md,q,d,w,U0,V,E,xv,hbar,LE,LC,E_0)

% Potential (global Ux och lokal U_v, mx m_v pss)
    [Ux,U_v,mx_tot,m_v, x_vloc,NL,NR, sb1, sb2, sw,d1start,wstart,d2start,d2fin] = simpUx_mx(xv,U0,V,0,w,d,1,1,1,1,md,mw);
    
    Ux = Ux+E_0;
    U_v = U_v+E_0;
%     E = E+E_0;

%     [Ux, U_v, mx_tot, m_v, x_loc, NL, NR,IDd1,IDw,IDd2] = ... 
%                             simpUx_mxV5(xv,U0,V,0,w,d,LE,LC,q,md,mw,E_0);
%     d1start = IDd1(1);
%     wstart = IDd1(2)+1;
%     d2start = IDd2(1);
%     d2fin = IDd2(2);
    
    % Övergång mellan emitter och dubbelbarriär
    xstartb1 = xv(d1start);
    xendb1 = xv(wstart-1);
    V_at_xstartb1 = Ux(d1start);
    V_at_xendb1 = Ux(wstart-1);
    %beta3 = -sb1;
    %alpha3 = V_at_xstartb1-beta3*xstartb1;
    emitterV = U_v(1);
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
    xstartw = xv(wstart);
    xendw = xv(d2start-1);
    V_at_xstartw = Ux(wstart);
    V_at_xendw = Ux(d2start-1);
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
    xstartd2 = xv(d2start);
    xendd2 = xv(d2fin);
    V_at_xstartd2 = Ux(d2start);
    V_at_xendd2 = Ux(d2fin);
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
    collectorV = U_v(end);
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

    
end

