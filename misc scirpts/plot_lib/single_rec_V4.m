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

x0 = -100;
xfin = 100;

Nx = 5000;
dx = (xfin-x0) / (Nx);
xv = (x0:dx:xfin)';
L = xfin-x0;

t0 = 0;
tfin = 9;
dt = 0.05e-2;
tv = (t0:dt:tfin)';
Nt = length(tv);

stab = dt/(dx^2)



xbar = -0;                      % Mitten av vågfunk. vid t=0
k0 = 2;                        % Medel vågtal/energi vid t=0, kopplad till r.mängd via p=hbar*k
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

nu = 100;
% Ux = -hbar^2/(2*m)*nu*(nu-1)*sech(xv-2);
% Ux = zeros(Nx+1,1);
Ux = 1.1*E*(heaviside(xv-15) - heaviside(xv-20));

figure(1)
plot(xv, Ux)
hold on
plot(xv, abs(psi0).^2+E)
hold off

%%

RK4 = Rk4_solve(psi0, Ux, hbar, m, k0, dx, dt, Nt);
% CN = CN_solve(psi0, Ux, Nx, Nt, dx, dt, m);
%%
nf = 4*60;

figure("Position",[200,200,1100,700])
skip = floor((Nt-1)/nf);
for j = 1:nf

    Tj = sprintf('%.5f',trapz(xv(3000:end), abs(RK4(3000:end,skip*j)).^2 ));
    tj = sprintf('%.4f',skip*j*dt);
    yyaxis left
    plot(xv, abs(RK4(:,skip*j)).^2)
    ylabel('$|\Psi|^2 \ [\mathrm{eV}/(\hbar c)]$','fontsize',14,'Interpreter','latex')
    ylim([-0,0.5])
%     hold on
    yyaxis right
    plot(xv, Ux)
    yline(E,'--',Color="#D95319")
    ylabel('$U(x), \ \langle E \rangle \ [\mathrm{eV}]$','fontsize',14,'Interpreter','latex')
    xlabel('$x \ [\mathrm{eV} \cdot \hbar c]$','fontsize',14,'Interpreter','latex')
%     plot(xv, abs(CN(:,skip*j)).^2, 'k-')
%     plot(xv, real(Psi1_ana(:,j)))
%     hold off
    legend('$|\Psi|^2$','$U(x)$', '$\langle E \rangle$','fontsize',12,'Interpreter','latex')
    ylim([-0,5])
%     xlim([-10,35])
    title(['T=', num2str(Tj), '        t=', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)

    
    drawnow
%     frame = j
    F(j)=getframe(gcf);
end


%%
nf = 600;

figure(3, "Position",[10,10,500,200])
skip = floor((Nt-1)/nf);
for j = 1:nf
    
    plot(p,abs( fftshift(fft(RK4(:,skip*j))) ).^2)
    hold on
%     plot(p,abs( fftshift(fft(RK4(:,1))) ),'--')
%     plot(xv, real(Psi1_ana(:,j)))
    hold off

    ylim([-2,60])
    xlim([-20,20])

    
    drawnow
%     frame = j
    
end


%%
video=VideoWriter('sgl_rec_v4');
video.Quality=100;
video.FrameRate=60;
open(video);
writeVideo(video,F);
close(video);