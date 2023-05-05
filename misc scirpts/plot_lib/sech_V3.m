close all
clear
clc
%%
% V = 1e-6; Evec = linspace(0,3*V,1000); w = 7/(sqrt(2*m*V));
c = 3e8;
eV = 1.602e-19;
hbarsc = 6.63e-34 / (2*pi);
m_e = 9.12e-31;

hbar = 1;
m = 1/2;

x0 = -30;
xfin = 70;

Nx = 2800;
dx = (xfin-x0) / (Nx);
xv = (x0:dx:xfin)';
L = xfin-x0;

t0 = 0;
tfin = 4.6;
dt = 0.01e-2;
tv = (t0:dt:tfin)';
Nt = length(tv);

stab = dt/(dx^2)



xbar = -10;                      % Mitten av vågfunk. vid t=0
k0 = 4;                        % Medel vågtal/energi vid t=0, kopplad till r.mängd via p=hbar*k
E = ((hbar*k0)^2) / (2*m);      % k översatt till kinetisk energi se 2.4 i Griffiths
sigmax = 10e-1;                   % Standardavvikelse i x vid t=0
sigmak = 1/(2*sigmax);          % Standardavvikelse i k vid t=0
aa = 1/(4*sigmax^2);            % Konstant från uppgift 2.21 Griffiths
vg = hbar*k0/m;                 % Grupphastighet Verkar stämma!

T_ana = @(E,V,w) (1 + (V.*sinh((w/hbar).*sqrt(2.*m.*(V-E))) ).^2./(4.*E.*(V-E))).^(-1);

dp = 2*pi/(dx*(Nx+1));
p = (-Nx/2:Nx/2)*dp;

% Normaliserat vågpaket som initialvillkor för numerisk lösare
Anorm = 1/(2*pi*sigmax^2)^(1/4);
psi0 = Anorm*exp(-aa*(xv-xbar).^2) .* exp(1i*k0*(xv));
normtest = trapz(xv, abs(psi0).^2);

nu = 20;
Ux = -hbar^2/(2*m)*nu*(nu-1)*sech(xv-2);
% Ux = zeros(Nx+1,1);
% Ux = 1.1*E*(heaviside(xv-15) - heaviside(xv-20));

figure(1)
plot(xv, Ux)
hold on
plot(xv, abs(psi0).^2)
hold off

%%

RK44 = Rk4_solve(psi0, Ux, hbar, m, k0, dx, dt, Nt);
% RK4 = CN_solve(psi0, Ux, Nx, Nt, dx, dt, m);
%%
nf = 500;

figure("Position",[200,200,1100,700])
skip = floor((Nt-1)/nf);
for j = 1:nf

    Tj = sprintf('%.5f',trapz(xv(1176:end), abs(RK44(1176:end,skip*j)).^2 ));
    tj = sprintf('%.4f',skip*j*dt);
    yyaxis left
    plot(xv, abs(RK44(:,skip*j)).^2)
    ylabel('$|\Psi|^2 \ [\mathrm{eV}/(\hbar c)]$','fontsize',14,'Interpreter','latex')
    ylim([-0.5,0.5])
%     hold on
    yyaxis right
    plot(xv, Ux)
    yline(E,'--',Color="#D95319")
    ylabel('$U(x), \ \langle E \rangle \ [\mathrm{eV}]$','fontsize',14,'Interpreter','latex')
    xlabel('$x \ [\mathrm{eV} \cdot \hbar c]$','fontsize',14,'Interpreter','latex')
%     plot(xv, abs(CN(:,skip*j)).^2, 'k-')
%     plot(xv, real(Psi1_ana(:,j)))
%     hold off
    legend('$|\Psi|^2$','$U(x)_{\mathrm{sech}}, \ \nu=20$', '$\langle E \rangle$','fontsize',12,'Interpreter','latex')
    ylim([-17,17])
    xlim([-20,20])
    title(['T=', num2str(Tj),',', '        t=', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)

    
    drawnow
%     frame = j
    F(j)=getframe(gcf);
end


%%

nf = 500;

figure("Position",[200,200,1100,700])
skip = floor((Nt-1)/nf);
for j = 1:nf

    Tj = sprintf('%.5f',trapz(xv(1176:end), abs(RK4(1176:end,skip*j)).^2 ));
    tj = sprintf('%.4f',skip*j*dt);
    yyaxis left
    plot(xv, real(RK44(:,skip*j)))
    hold on
    plot(xv, imag(RK44(:,skip*j)))
    hold off
    ylabel('$\Psi (x,t)$','fontsize',14,'Interpreter','latex')
    ylim([-1,1])

    yyaxis right
    plot(xv, Ux)
 
    ylabel('$U(x) \ [\mathrm{eV}]$','fontsize',14,'Interpreter','latex')
    xlabel('$x \ [\mathrm{eV} \cdot \hbar c]$','fontsize',14,'Interpreter','latex')

    legend('$\mathrm{Re}(\Psi)$', '$\mathrm{Im}(\Psi)$','$U(x)_{\mathrm{sech}}, \ \nu=20$','fontsize',12,'Interpreter','latex')
    ylim([-17,17])
    xlim([-20,20])
    title(['T=', num2str(Tj),',', '        t=', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)

    
    drawnow
%     frame = j
    F(j)=getframe(gcf);
end

%%
video=VideoWriter('sech_V3_ReIm');
video.Quality=100;
video.FrameRate=60;
open(video);
writeVideo(video,F);
close(video);

%% Se boken "Spectral Methods in Matlab" i onedrive>KEX_delat>lit>pseudo-spec
nline = 24;
skip_line = floor((Nt-1)/nline)
%%
figure("Position",[200,200,1100,700])
waterfall(xv,tv(1:skip_line:Nt-1),abs(RK44(:,1:skip_line:Nt-1).').^2)
view(10,70) 
% colormap([0 0 0])
axis([x0 xfin 0 tfin 0 0.4]) 
ylabel('$t =  [\hbar/ \mathrm{eV}]$','fontsize',10,'Interpreter','latex')
zlabel |\Psi|^2, grid on, xlabel x
print('sech_V3_wf_c','-depsc')
print('sech_V3_wf_c','-dpng')
%%
nlinesurf = 400;
skip_linesurf = floor((Nt-1)/nlinesurf)
figure("Position",[200,200,1100,700])
surface(xv(1:5:Nx),tv(1:skip_linesurf:Nt-1),abs(RK44((1:5:Nx),1:skip_linesurf:Nt-1).').^2, 'EdgeColor','none')
cb = colorbar;
hColourbar.Label.Position(1) = 3;
ylabel(cb, '$|\Psi|^2$','fontsize',14,'Interpreter','latex')
ylabel('$t =  [\hbar/ \mathrm{eV}]$','fontsize',14,'Interpreter','latex')
xlabel('$x \ [\mathrm{eV} \cdot \hbar c]$','fontsize',14,'Interpreter','latex')
title('$U(x)_{\mathrm{sech}}, \ \nu=20$', '$\langle E \rangle=16 \ [\mathrm{eV}]$','fontsize',14,'Interpreter','latex')
xlim([x0,xfin])
ylim([0,tfin])
print('sech_V3_surf','-depsc')
print('sech_V3_surf','-dpng')