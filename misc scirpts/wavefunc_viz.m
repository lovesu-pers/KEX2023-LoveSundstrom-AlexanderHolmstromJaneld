clear 
clc
close all
set(0, 'DefaultLineLineWidth', 3);
C = linspecer(4);
RTD_input;


Nx = 8000;%5193;
NVb = 700; 
Nt = 1e3;
Nframe = 600;

x0 = -500;
xfin = L+500;
dx = (xfin-x0) / (Nx+1);
d = dx*round(d/dx);
w = dx*round(w/dx);
% dxTMM = (L) / (NxTMM+1)
xv = linspace(x0,xfin,Nx+2);

V0 = 0;
Vfin = 0.5/(V_Hartree_SI);
Vbv = linspace(V0,Vfin,NVb);


Ux_mat = zeros(Nx+2,NVb);
for Vj = 1:NVb
    [Ux_mat(:,Vj),~,mx ,~,~ ,~,NR,~,~,~] = ... 
                            simpUx_mxV5(xv,(U0),Vbv(Vj),xb,w,d,LE,LC,q,md,mw,E_0);
end

xbar = -65*nm_a0_sf;
EvCN = linspace(E_0+beta*1e-5,EvE(end),NE);



kvec = sqrt(2*mw*(EvCN) )/(hbar);
k0 = kvec(1);
kfin = kvec(end);
S = L+abs(2*xbar);
tmax = S/(hbar*k0/mw)
tmin = S/(hbar*kfin/mw)
int_index = (Nx-NR);


sigmax = 5*d;
aa = 1/(4*sigmax^2); 

% Anorm = 1/(2*pi*sigmax^2)^(1/4);
% psi0 = Anorm*exp(-aa*(xv(2:end-1).'-xbar).^2) .* exp(1i*k0*(xv(2:end-1).'));
% normtest = trapz(xv(2:end-1), abs(psi0).^2);
% kspace = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
% dk = 2*pi/(dx*Nx);
% kv = (-Nx/2 : Nx/2-1)*dk;

figure(1)
plot(xv/nm_a0_sf,Ux_mat(:,1)*Eh_eV)
hold on
plot( xv/nm_a0_sf,Ux_mat(:,end)*Eh_eV,':')
% plot(xv(2:end-1)/nm_a0_sf,abs(psi0)+EvE(1)*Eh_eV)
hold off
yline(E_FE*Eh_eV,'--')
xlabel('x [nm]')
ylabel('Potential profile [eV]')
% legend('V_b = 0 [V]', ['V_b = ' num2str(Vfin*V_Hartree_SI) ' [V]'], 'E_F^E')
xlim([-50,50])

% figure(2)
% plot(kv, abs(kspace))

TvQ_EC = zeros(NVb,1);
TvQ_CE = zeros(NVb,NE);
PSI_EC = zeros(Nx+2,NVb);
PSI_CE = zeros(Nx+2,NE);
cx = zeros(Nx+2,NVb);
jvQ = zeros(NVb,1);
%%
for jQ = 1:NVb
    VjQ = Vbv(jQ);
    E_FC = E_FE-VjQ*q;

    EvC = linspace(E_0 - VjQ*q - 1e-6, E_FC+20*kB*T, NE);
    
    
    EkQ = 0.0540;
    [PSI_EC(:,jQ), TvQ_EC(jQ)] = QTBM_func_Hartree_V2(xv,U0,VjQ, ... 
        EkQ,xb,w,d,LE,LC,md,mw,E_0,1,0,0);
    
    
    
       

    

    
%     eedens = TvQ_EC(jQ,:).*supply_function_fd(EvE,E_FE,E_FC,kB,T);
  

%     jvQ(jQ) = (mw*q*beta/(2*pi^2*hbar^3)) * trapz(EvE,eedens);
    
%     fE = log(1+exp( (E_FE-EvE.') / (beta)));
%     fC = log(1+exp( (E_FC-EvC.') / (beta)));

%     cEC = trapz(EvE,abs(PSI_EC).^2 .* fE.',2);
%     cCE = trapz(EvC,abs(PSI_CE).^2 .* fC.',2);
%     
%     cx(:,jQ) = (mw*q*beta/(2*pi^2*hbar)) * (cEC+cCE);

    disp(['QTBM bias point: ' num2str(jQ) '/' num2str(NVb) ]) 
end
%%
set(groot, 'DefaultAxesFontSize',16)
set(groot,"defaultFigurePosition", [10,10,900,800])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
C = linspecer(3);

[~, V2_index] = min(abs(Vbv*V_Hartree_SI-0.01))
[~, V3_index] = min(abs(Vbv*V_Hartree_SI-0.048641))
[~, V4_index] = min(abs(Vbv*V_Hartree_SI-0.125))




f4 = figure(Position=[10,10,900,800]);
yyaxis left
plot(xv/nm_a0_sf,Ux_mat(:,1)*Eh_eV,'-.')
yline(EkQ*Eh_eV,':',LineWidth=3)
ylim([1.2,1.65])
ylabel('Potential profile [eV]')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
ax1 = gca;
ax1.TickLabelInterpreter = "latex";
ax1.FontSize = 18;


yyaxis right
plot(xv/nm_a0_sf, abs(PSI_EC(:,1)).^2)
ylim([-4.5,4.5])
ylabel('$|\psi(x)|^2$ [1/nm]')
xlim([-15,20])
ax2 = gca;
grid on
title('$V_b=0$ [mV]')
legend('$V(x)$','$E$','$|\psi(x)|^2$')
print(f4,'psiQTBM1', '-depsc')
%%
f5 = figure(Position=[10,10,900,800]);
yyaxis left
plot(xv/nm_a0_sf,Ux_mat(:,V2_index)*Eh_eV,'-.')
yline(EkQ*Eh_eV,':',LineWidth=3)
ylim([1.2,1.65])
ylabel('Potential profile [eV]')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
ax1 = gca;
ax1.TickLabelInterpreter = "latex";
ax1.FontSize = 18;


yyaxis right
plot(xv/nm_a0_sf, abs(PSI_EC(:,V2_index)).^2)
ylim([-4.5,4.5])
ylabel('$|\psi(x)|^2$ [1/nm]')
xlim([-15,20])
ax2 = gca;
grid on
legend('$V(x)$','$E$','$|\psi(x)|^2$')
title('$V_b=10$ [mV]')
print(f5,'psiQTBM2', '-depsc')

%%
f6 = figure(Position=[10,10,900,800]);
yyaxis left
plot(xv/nm_a0_sf,Ux_mat(:,V3_index)*Eh_eV,'-.')
yline(EkQ*Eh_eV,':',LineWidth=3)
ylim([1.2,1.65])
ylabel('Potential profile [eV]')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
ax1 = gca;
ax1.TickLabelInterpreter = "latex";
ax1.FontSize = 18;


yyaxis right
plot(xv/nm_a0_sf, abs(PSI_EC(:,V3_index)).^2)
ylim([-4.5,4.5])
ylabel('$|\psi(x)|^2$ [1/nm]')
xlim([-15,20])
ax2 = gca;
grid on
legend('$V(x)$','$E$','$|\psi(x)|^2$')
title('$V_b=48$ [mV]')
print(f6,'psiQTBM3', '-depsc')

%%

f7 = figure(Position=[10,10,900,800]);
yyaxis left
plot(xv/nm_a0_sf,Ux_mat(:,V4_index)*Eh_eV,'-.')
yline(EkQ*Eh_eV,':',LineWidth=3)
ylim([1.2,1.65])
ylabel('Potential profile [eV]')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
ax1 = gca;
ax1.TickLabelInterpreter = "latex";
ax1.FontSize = 18;


yyaxis right
plot(xv/nm_a0_sf, abs(PSI_EC(:,V4_index)).^2)
ylim([-4.5,4.5])
ylabel('$|\psi(x)|^2$ [1/nm]')
xlim([-15,20])
ax2 = gca;
grid on
legend('$V(x)$','$E$','$|\psi(x)|^2$')
title('$V_b=125$ [mV]')
print(f7,'psiQTBM4', '-depsc')
%%
figure(3)
plot(TvQ_EC)
%%

figure(4)
plot(Vbv*V_Hartree_SI,TvQ_EC)
%%
psimax = max(max(abs(PSI_EC)));
figure(5)
yyaxis left
psifig = plot(xv/nm_a0_sf, real(PSI_EC(:,1)),LineWidth=1);
% ylim([0,7e-5])

yyaxis right
plot(xv/nm_a0_sf,1*Ux_mat(:,1)+0.01,'--',LineWidth=1)
grid on

xlim([-15,20])
%%
for fj = 1:NVb
    psifig.YData = abs(PSI_EC(:,fj));
    pause(0.001)
    drawnow
end
