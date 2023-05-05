close all
clear
clc

%%
hbar = 1;
m = 1;

x0 = -30;
xfin = 30;

Nx = 1000;
dx = (xfin-x0) / (Nx);
xv = (x0:dx:xfin)';
L = xfin-x0;

t0 = 0;
tfin = 1;
dt = 0.25e-3;
tv = (t0:dt:tfin)';
Nt = length(tv);

stab = dt/(dx^2);


kappa = 0.1; %1/40;

xbar = pi;                      % Mitten av vågfunk. vid t=0
E = 0;
k0 =sqrt( E*2*m )/ hbar;                        % Medel vågtal/energi vid t=0, kopplad till r.mängd via p=hbar*k

sigmax = 0.25;                   % Standardavvikelse i x vid t=0
sigmak = 1/(2*sigmax);          % Standardavvikelse i k vid t=0

aa = 1/(4*sigmax^2);            % Konstant från uppgift 2.21 Griffiths
vg = hbar*k0/m;                 % Grupphastighet Verkar stämma!



% Normaliserat vågpaket som initialvillkor för numerisk lösare

Anorm = 1/(2*pi*sigmax^2)^(1/4);



% figure(2)
% plot(kv, abs(k0_v).^2)


% kprob = trapz(kv, kv.*abs(k0_v).^2)

lambda = 6;
De = 20;
alpha = -0.16;
g = 100;
% Ux = 2*(0.1*xv.^6 + 4*xv.^3 - 5*xv.^2 - 4*xv);
% Ux = -0.5*lambda*(lambda+1)*(sech(xv)).^2;
% Ux = 10*E*(heaviside(xv+10) - heaviside(xv+9)) + 10*E*(heaviside(xv-39) - heaviside(xv-40));
% Ux = 0.5*kappa*xv.^2;
% Ux = 0.25*kappa*xv.^4 - 5*kappa*xv.^2;
% Ux = kappa*xv.^4;
Ux =  m*g*(1-cos(xv));
Ux_pend = m*g*(xv).^2;
% Ux = De * (1-exp(alpha*xv)).^2;

H = @(x,k) (k.^2)/2 + m*g*(1-cos(x));

xtemp = -10:0.01:10;
ktemp = -20:0.1:20;
[xtest, ktest] = meshgrid(xtemp, ktemp);

Htest = H(xtest, ktest);

% figure(1)
% plot(xv, abs(psi0).^2, xv, Ux, '--')

[egenfuncs, egenv, conv] = egenfunc(Nx+1, Ux_pend, m, hbar, dx, 8);

conv
% psi0 = egenfuncs(:,2);
psi0 = Anorm*exp(-aa*(xv-xbar).^2) .* exp(1i*k0*(xv));

figure(2)
plot(xv, 50*egenfuncs(:,4))
hold on
plot(xv, 10*abs(psi0).^2, xv, Ux, '--',xv,Ux_pend)
yline(egenv(1,1))
yline(egenv(2,2))
yline(egenv(3,3))
hold off
ylim([-50, 100])
xlim([-6,6])



normtest = trapz(xv, abs(psi0).^2);


k0_v = fftshift( (1/sqrt(2*pi))*fft(psi0(1:Nx))*dx );
dk = 2*pi/(dx*Nx);
kv = ((-Nx/2:Nx/2-1)*dk).';
Nk = length(kv);

figure(3)
plot(kv, abs(k0_v).^2)

figure(4)
contour(xtemp,ktemp, Htest,0:20:300)

%%
CN = CN_solve(psi0,Ux,Nx,Nt,dx,dt,m);
%%
close all
nf = 100;
skip_frame = floor((Nt-1) / nf);
figure('Position',[100, 100, 1000, 900])
F = struct('cdata', cell(1,nf), 'colormap', cell(1,nf));
for j = 1:nf
    tj = sprintf('%.4f',skip_frame*j*dt);
    subplot(2,1,1)

    plot(xv, Ux, '-k')
    hold on
    plot(xv, 10*abs( CN(:,skip_frame*j) ).^2  )
    xlabel('$x \ [\mathrm{eV}^{-1} \cdot \hbar c]$','fontsize',15,'Interpreter','latex')
    legend('$U(x)$','$|\Psi(t)|^2$','fontsize',15,'Interpreter','latex')
    title('$ U(x)=mg(1-\mathrm{cos}(x)), \ m=1 \ g=100$','fontsize',18,'Interpreter','latex')
    subtitle(['        t=', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)
    ylim([0, 50])
    xlim([-7,7])
    hold off
    
    subplot(2,1,2)
    plot(kv(1:Nx), abs(fftshift( (1/sqrt(2*pi))*fft(CN(1:Nx, j*skip_frame))*dx )))
    xlabel('$k \ [\mathrm{eV} \cdot \hbar^{-1}c^{-1} ]$','fontsize',15,'Interpreter','latex')
    ylabel('$|\Phi(k)|^2$','fontsize',15,'Interpreter','latex')
    xlim([-50,50])
    ylim([-0,1])
    
    drawnow
    F(j) = getframe(gcf);
end

%%
tic
Wv = zeros(Nx,Nk,nf);
for jj = 1:nf
    Wj = mywigner(CN(1:Nx,skip_frame*jj));
    Wv(:,:,jj) = Wj ;
end
toc
W_U = mywigner(Ux(1:Nx));
%%
close all
nf = 100;
nl = 500;
% skip_l = floor((Nt-1) / Nk);
skip_l = 1;
skip_frame = round((Nt-1) / nf);

figure('Position',[0, -100, 1000, 900])
F2 = struct('cdata', cell(1,nf), 'colormap', cell(1,nf));

zmax = max(max(max(Wv)));
zmin = min(min(min(Wv)));
for j = 1:nf
    tj = sprintf('%.4f',skip_frame*j*dt);
    Wj = Wv(1:skip_l:end,1:skip_l:end,j).';
    
    contour3( xtemp, ktemp, Htest,0:20:300,'-b'   )
    hold on
    surf(xv(1:skip_l:Nx), kv(1:skip_l:Nx), Wj, 'EdgeColor','none', 'FaceColor','interp')
    hold off
    
    xlim([-2,8])
    ylim([-20,20])
%     zlim([zmin,20*zmax])
    zlim([-30,420])

    ylabel('$k \ [\mathrm{eV} \cdot \hbar^{-1}c^{-1} ]$','fontsize',15,'Interpreter','latex')
    xlabel('$x \ [\mathrm{eV}^{-1} \cdot \hbar c]$','fontsize',15,'Interpreter','latex')

    title('$\mathrm{Wigner-function}, \ U(x)=mg(1-\mathrm{cos}(x)), \ m=1 \ g=100$','fontsize',18,'Interpreter','latex')
    subtitle(['        t=', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)

    view(2)
    %     view(-40, 45) 

    colormap(flipud(cbrewer2('RdGy')));
    clim([-10, 30])
    
    c = colorbar;
    ylabel(c,'$W(x,k)$','fontsize',15,'Interpreter','latex');
    
    drawnow
    F2(j) = getframe(gcf);
end

%%
close all

nn = 5;
mm = 1;

Nt = 500;
tv = linspace(0,3,Nt);
figure(3)

psin = egenfuncs(:,nn);
En = egenv(nn,nn);

psim = egenfuncs(:,mm);
Em = egenv(mm,mm);

sup = zeros(Nx,Nt);

for tt = 1:Nt
    mn = psin*exp(-En*1i*tv(tt)/hbar)+ psim*exp(-Em*1i*tv(tt)/hbar);
    sup(:,tt) = mn;
    plot(xv(1:Nx), real(mn))
    hold on
    plot(xv(1:Nx), abs(mn))
    plot(xv(1:Nx), imag(mn))
    hold off
    ylim([-0.3,0.3])
    drawnow
end
%%
video=VideoWriter('QP_Wigner_2D_V6','MPEG-4');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,F2);
close(video);

%%
video=VideoWriter('QP_psi_phi_V6','MPEG-4');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,F);
close(video);
