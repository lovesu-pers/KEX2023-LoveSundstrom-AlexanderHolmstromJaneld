clear 
clc
set(0, 'DefaultLineLineWidth', 1);
%%
RTD_input;

md = mw;

d = 20;
w = 0;
x0 = -2000;
xfin = 2*d+3000;


% Nxv = [447;644;742;841;1432;1727;2121;2614;2811;3500;4386;5962;9114];
% Nxv = [447; 1727; 2614; 5962; round(linspace(6000,8000,20)')];
% Nxv = [250;300;350;400;447;550;600;644;700;742;800;841;900;1000;1432;1500;1727;1800;2121;2500;2614;2700;2811;3000;3500;4000;4386;4500;5000;5500;5962;9114];
Nxv = [447 550 700 800 1432 1800 2811 5500].';
% dxv = flip([0.02 0.08 0.1 0.2 0.25 0.5 1 2 4 8 16]);
% Nxv = round((xfin-x0)./(dxv)-1);
% Nxv = round((linspace(1000,15000,100)).');
NN = length(Nxv);
Nx = Nxv(1);
NE = 1;
E = 0.05/Eh_eV;


NVb = 1; 
Nt = 1e3;
Nframe = 600;


dx = (xfin-x0) / (Nx+1);


% dxTMM = (L) / (NxTMM+1)
xv = linspace(x0,xfin,Nx+2);

V0 = 0;
Vfin = 0.5/(V_Hartree_SI);
Vbv = 0;



    [Ux,~,mx ,~,~ ,~,NR,~,~,~] = ... 
                            simpUx_mxV5(xv,(U0),Vbv,xb,w,d,LE,LC,q,md,mw,E_0);

EvCN = linspace(E_FE*0.5,U0,NE);

xbar = -45*nm_a0_sf;

k0 = sqrt(2*mw*(E-E_0) )/(hbar);


S = L+abs(2.5*xbar);
tmax = S/(hbar*k0/mw);
int_index = (Nx-NR)+1;

[1,2,3,4,5];
sigmax = 1*d;
aa = 1/(4*sigmax^2); 

Anorm = 1/(2*pi*sigmax^2)^(1/4);
psi0 = Anorm*exp(-aa*(xv(2:end-1).'-xbar).^2) .* exp(1i*k0*(xv(2:end-1).'));
normtest = trapz(xv(2:end-1), abs(psi0).^2);
kspace = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
dk = 2*pi/(dx*Nx);
kv = (-Nx/2 : Nx/2-1)*dk;

figure(1)
plot(xv/nm_a0_sf,Ux*Eh_eV)
hold on
plot( xv/nm_a0_sf,Ux*Eh_eV,':')
plot(xv(2:end-1)/nm_a0_sf,abs(psi0)+EvE(1)*Eh_eV)
hold off
yline(E_FE*Eh_eV,'--')
xlabel('x [nm]')
ylabel('Potentail profile [eV]')
xlim([-50,50])

figure(2)
plot(kv, abs(kspace))

% 
% aerrCN = zeros(NN,1);
% 
% TvQ = zeros(NN,1);
% TvAna = T_Tong(E,U0,2*d,mw,hbar);
% TvTMMC = zeros(NN,1);
% TvTMML = zeros(NN,1);
% TvCN = zeros(NN,1);

%%
TvCN5 = zeros(NN,1);
for Nj = flip(1:NN)

    disp(['loop:' num2str(Nj), '/', num2str(NN)])
    Nx = Nxv(Nj);
    dx = (xfin-x0) / (Nx+1);
    xv = linspace(x0,xfin,Nx+2);

    [Ux,~,mx ,~,~ ,~,NR,~,~,~] = ... 
                            simpUx_mxV5(xv,(U0),Vbv,xb,w,d,LE,LC,q,md,mw,E_0);
    int_index = (Nx-NR)+1;




    tkCN = S/(hbar*k0/mw);
    dt = tkCN/(Nt);
    tv = linspace(0,tkCN,Nt);

    psi_kCN = Anorm*exp(-aa*(xv(2:end-1).'-xbar).^2) .* exp(1i*k0*(xv(2:end-1).'));
    tic
    PSI_CN = CN_solve_mx(psi_kCN, Ux(2:end-1), mx(2:end-1), Nx, Nt, Nframe, dx, dt, hbar,0);
    TvCN5(Nj) = trapz((xv(int_index:Nx)).', abs(PSI_CN(int_index:end)).^2);
    tj = toc;
    disp(['Last time: ', num2str(tj)])


end
%%
set(groot, 'DefaultAxesFontSize',16)
set(groot,"defaultFigurePosition", [10,10,760,590])
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

dxv = (xfin-x0) ./ (Nxv+1);
T_ana = T_Tong(E,U0,2*d,mw,hbar);
err = abs(TvCN-T_ana);
C = linspecer(3);

figure(2)
loglog(dxv, err,'o')
legend('CN', ...
    'Location','northwest','NumColumns',2)
xlabel('Grid size, $\Delta x$ [nm]','Interpreter','latex')
ylabel('Absolute error, $|T_{\mathrm{analytic}}-T_{\mathrm{CN}}|$ [-]')
grid on
% print(2,'conv_CN2', '-depsc')


figure(3)
loglog(Nxv,err,'*')




%%
psimax = max(max(abs(PSI_CN)));
figure(5)
psifig = plot(xv(1:Nx)/nm_a0_sf, abs(PSI_CN(:,1)),LineWidth=1);
hold on
plot(xv/nm_a0_sf,1*Ux+0.01,'--',LineWidth=1)
hold off
grid on
ylim([-0,0.024])
xlim([-50,50])
%%
for fj = 1:Nframe
    psifig.YData = abs(PSI_CN(:,fj));
    pause(0.01)
    drawnow
end

%%
T_ana = T_Tong(E,U0,2*d,mw,hbar);
err1 = abs(TvCN1-T_ana);
err2 = abs(TvCN2-T_ana);
err3 = abs(TvCN3-T_ana);
err4 = abs(TvCN4-T_ana);
err5 = abs(TvCN5-T_ana);

figure(6)
loglog(Nxv,err1)
hold on
% loglog(Nxv,err2)
% loglog(Nxv,err3)
% loglog(Nxv,err4)
loglog(Nxv,err5)
grid on
hold off
legend('$\delta_x=5d$','$\delta_x=d$')
xlabel('Space steps, $N_x$ [-]','Interpreter','latex')
ylabel('Absolute error, $|T_{\mathrm{analytic}}-T_{\mathrm{CN}}|$ [-]')
% print(6,'conv_deltax_CN', '-depsc')

%%

Nx = 2811;
xv = linspace(x0,xfin,2811+2);
dx = xv(2)-xv(1);


sigmax = 5*d;
aa = 1/(4*sigmax^2); 
Anorm = 1/(2*pi*sigmax^2)^(1/4);
psi0 = Anorm*exp(-aa*(xv(2:end-1).'-xbar).^2) .* exp(1i*k0*(xv(2:end-1).'));
normtest = trapz(xv(2:end-1), abs(psi0).^2);
kspace1 = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
dk = 2*pi/(dx*Nx);
kv = (-Nx/2 : Nx/2-1)*dk;

sigmax = 1*d;
aa = 1/(4*sigmax^2); 
Anorm = 1/(2*pi*sigmax^2)^(1/4);
psi0 = Anorm*exp(-aa*(xv(2:end-1).'-xbar).^2) .* exp(1i*k0*(xv(2:end-1).'));
normtest = trapz(xv(2:end-1), abs(psi0).^2);
kspace2 = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
dk = 2*pi/(dx*Nx);
kv = (-Nx/2 : Nx/2-1)*dk;

figure(7)
plot(kv/nm_a0_sf, abs(kspace1).^2*nm_a0_sf)
hold on
plot(kv/nm_a0_sf, abs(kspace2).^2*nm_a0_sf,'--')
hold off
legend('$\delta_x=5d$','$\delta_x=d$')
xlabel('Wave number, $k$ [1/nm]','Interpreter','latex')
ylabel('Wave function, $|\phi(k)|^2$ [nm]')
xlim([-4e-3,4e-3])
% print(7,'kspace_deltax', '-depsc')

trapz(kv/nm_a0_sf,abs(kspace2).^2*nm_a0_sf)