clear

RTD_input;
md = md;


% w = 2.35*w 
w = 1*w;
d = 1.5*d;
Nx = 2500;
NVb = 25; 
Nt = 10e3;
Nframe = 800;

Ec = E_0;
% Ec = 0;

x0 = -5000;
xfin = 5000;

dx = (xfin-x0) / (Nx+1);

xv = (linspace(x0,xfin,Nx))';

Vb = 0%0.150/(V_Hartree_SI);
% Vb = 0.0030624;




[Ux,~,mx ,~,~ ,~,NR,~,~,~] = ... 
                        simpUx_mxV5(xv,(1*U0),Vb,xb,w,d,LE,LC,q,md,mw,Ec);

xbar = -1000;

%E0 = 0.0125; %E_FE-Ec;
% E0 = 0.00343671;
E0 = 0.00308236;
%


k0 =  sqrt(2*mw*(E0) )/(hbar)

vg1 =  hbar*sqrt(2*mw*(E0) )/(hbar*mw);
vg2 =  hbar*sqrt(2*mw*(E0-Ec) )/(hbar*mw);



tmax = ceil((2000-xbar)/(vg1));

% dt = tmax/(Nt);
tv = linspace(0,tmax,Nt);
dt = tv(2)-tv(1)

% vg2*tmax

sigmax = 3*d;
aa = 1/(4*sigmax^2); 

Anorm = 1/(2*pi*sigmax^2)^(1/4);
psi0 = Anorm*exp(-aa*(xv-xbar).^2) .* exp(1i*k0*(xv));
normtest = trapz(xv, abs(psi0).^2);

kspace = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
dk = 2*pi/(dx*Nx);
kv = (-Nx/2 : Nx/2-1)*dk;

 (Nx/2-1)*dk

int_index = (Nx-NR)

figure(1)
yyaxis left
plot(Ux)
yyaxis right
plot(abs(psi0) )
ylim([0,0.50])
xlim([550 1000])


figure(2)
plot(kv,abs(kspace).^2)
%%

PSI_CN = CN_solve_mx(psi0, Ux, mx, Nx, Nt, Nframe, dx, dt, hbar,1);
%%
test = psi0;
PSIold = CN_solve(test, Ux, Nx-1, Nt, dx, dt, mw);

figure(70)
plot(xv,abs(PSIold(:,end)))
%%
set(0, 'DefaultLineLineWidth', 3);
C = linspecer(4);

psimax = max(max(abs(PSI_CN)));
figure(5)
yyaxis left
psifig = plot(xv/nm_a0_sf, abs(PSI_CN(:,10)),LineWidth=2);
hold on
realfig = plot(xv/nm_a0_sf, real(PSI_CN(:,10)),'-',LineWidth=2,Color=C(3,:));
hold off
ylim([-0.08,0.08])
yyaxis right
plot(xv/nm_a0_sf,1*Ux,'-',LineWidth=1)
xlim([-1000/nm_a0_sf,1000/nm_a0_sf])
grid on
% legend('$|\Psi|$','Re$(\Psi)$','$V(x)$', 'Location','northwest',fontsize=10)
xlabel('$x$ [nm]')


T = trapz(xv(int_index:end), abs(PSI_CN(int_index:end,end)).^2)
% ylim([-0,0.024])



%%
for fj = 1:770
    psifig.YData = abs(PSI_CN(:,fj));
    realfig.YData = real(PSI_CN(:,fj));
%     pause(0.01)
    fj
    drawnow
    F(fj) = getframe(gca);
end
%%
video=VideoWriter('qbs156t2_v1z');
video.Quality=100;
video.FrameRate=60;
open(video);
writeVideo(video,F);
close(video);
%%
figure(4)
waterfall(xv/nm_a0_sf,tv(1:20:770),abs(PSI_CN(:,1:20:770).'))
view(18,44) 
box off
axis off
colormap([0 0 0])
% xlim(axes1,[-299.334817678919 299.280488237295]);
% ylim(axes1,[0.273947583705669 798.427688805324]);
% zlim(axes1,[4.36593532509523e-05 0.0798590316914002]);
% axis([-150 150 0 180 0 0.1]) 
% ylabel('$t$ [au]'), zlabel $|\psi|$, grid off, xlabel('$x$ [nm]') 
grid off
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
set(gca,'ZTickLabel',[]);
set(gca,'Color','w')
box off

%%
figure(50)
waterfall(xv/nm_a0_sf,tv(1:20:770),abs(PSI_CN(:,1:20:770).'))
view(18,44) 
% colormap([0 0 0])
% axis([-150 150 0 200 0 0.1]) 
%ylabel('$t$ [au]'), zlabel $|\psi|$, grid on, xlabel('$x$ [nm]') 
%%
% close all
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot, 'DefaultAxesFontSize',16)
set(groot,"defaultFigurePosition", [10,10,900,800])


t1 = tmax/5*t_sf;
t2 = tmax*3/5*t_sf;
t3 = tmax*5/5*t_sf;


figure(6)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,1)))
ylabel('$|\psi(x)|$ [1/$\sqrt{\mathrm{nm}}$]')
xlim([-100,50])
ylim([-0.08,0.09])
yyaxis right
plot(xv/nm_a0_sf, Ux*Eh_eV,':')
yline(E_FE*Eh_eV,'--',LineWidth=1)
ylabel('Potential profile [eV]')
grid on
title('$t=0$ [s]')
legend('$|\psi(x)|$','$V(x)$','$E$')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
% print(6,'CN_wavefunc1', '-depsc')

figure(7)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,round(end/5))))
ylabel('$|\psi(x)|$ [1/$\sqrt{\mathrm{nm}}$]')
xlim([-50,50])
ylim([-0.08,0.09])
yyaxis right
plot(xv/nm_a0_sf, Ux*Eh_eV,':')
yline(E_FE*Eh_eV,'--',LineWidth=1)
ylabel('Potential profile [eV]')
grid on
title('$t=0.9\cdot 10^{-13}$ [s]')
legend('$|\psi(x)|$','$V(x)$','$E$')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
% print(7,'CN_wavefunc2', '-depsc')

figure(8)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,round(3*end/5))))
ylabel('$|\psi(x)|$ [1/$\sqrt{\mathrm{nm}}$]')
xlim([-10,20])
ylim([-0.02,0.05])
yyaxis right
plot(xv/nm_a0_sf, Ux*Eh_eV,':')
yline(E_FE*Eh_eV,'--',LineWidth=1)
ylabel('Potential profile [eV]')
grid on
title('$t=2.85\cdot 10^{-13}$ [s]')
legend('$|\psi(x)|$','$V(x)$','$E$')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
% print(8,'CN_wavefunc3', '-depsc')

figure(9)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,round(5*end/5))))
ylabel('$|\psi(x)|$ [1/$\sqrt{\mathrm{nm}}$]')
xlim([-10,20])
ylim([-0.028,0.05])
yyaxis right
plot(xv/nm_a0_sf, Ux*Eh_eV,':')
yline(E_FE*Eh_eV,'--',LineWidth=1)
ylabel('Potential profile [eV]')
grid on
title('$t=4.77\cdot 10^{-13}$ [s]')
legend('$|\psi(x)|$','$V(x)$','$E$')
xlabel('$x$ [nm]',Interpreter='latex',FontSize=16)
% print(9,'CN_wavefunc4', '-depsc')

%%
set(groot,"defaultFigurePosition", [10,10,500,500])
NN = 770;

t2 = round(NN*2/5);
t3 =  round(NN*0.55);
t4 =  round(NN*4/5);

figure(6)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,1)),LineWidth=2);
hold on
plot(xv/nm_a0_sf, real(PSI_CN(:,1)),'-',LineWidth=2,Color=C(3,:));
hold off
ylim([-0.1,0.1])
yyaxis right
plot(xv/nm_a0_sf,1*Ux,'-',LineWidth=1)
xlim([-1000/nm_a0_sf,1000/nm_a0_sf])
grid on
legend('$|\Psi|$','Re$(\Psi)$','$V(x)$', 'Location','northwest',fontsize=10)
xlabel('$x$ [nm]')
% print(6,'CN2_wavefunc1', '-depsc')

figure(7)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,t2)),LineWidth=2);
hold on
plot(xv/nm_a0_sf, real(PSI_CN(:,t2)),'-',LineWidth=2,Color=C(3,:));
hold off
ylim([-0.1,0.1])
yyaxis right
plot(xv/nm_a0_sf,1*Ux,'-',LineWidth=1)
xlim([-1000/nm_a0_sf,1000/nm_a0_sf])
grid on
legend('$|\Psi|$','Re$(\Psi)$','$V(x)$', 'Location','northwest',fontsize=10)
xlabel('$x$ [nm]')
% print(7,'CN2_wavefunc2', '-depsc')

figure(8)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,t3)),LineWidth=2);
hold on
plot(xv/nm_a0_sf, real(PSI_CN(:,t3)),'-',LineWidth=2,Color=C(3,:));
hold off
ylim([-0.1,0.1])
yyaxis right
plot(xv/nm_a0_sf,1*Ux,'-',LineWidth=1)
xlim([-1000/nm_a0_sf,1000/nm_a0_sf])
grid on
legend('$|\Psi|$','Re$(\Psi)$','$V(x)$', 'Location','northwest',fontsize=10)
xlabel('$x$ [nm]')
% print(8,'CN2_wavefunc3', '-depsc')

figure(9)
yyaxis left
plot(xv/nm_a0_sf, abs(PSI_CN(:,t4)),LineWidth=2);
hold on
plot(xv/nm_a0_sf, real(PSI_CN(:,t4)),'-',LineWidth=2,Color=C(3,:));
hold off
ylim([-0.1,0.1])
yyaxis right
plot(xv/nm_a0_sf,1*Ux,'-',LineWidth=1)
xlim([-1000/nm_a0_sf,1000/nm_a0_sf])
grid on
legend('$|\Psi|$','Re$(\Psi)$','$V(x)$', 'Location','northwest',fontsize=10)
xlabel('$x$ [nm]')
% print(9,'CN2_wavefunc4', '-depsc')
