close all
clear
clc

%%
hbar = 1;
m = 1;

x0 = -20;
xfin = 20;

Nx =800;
dx = (xfin-x0) / (Nx);
xv = (x0:dx:xfin)';
L = xfin-x0;

t0 = 0;
tfin = 2.2;
dt = 0.2e-3;
tv = (t0:dt:tfin)';
Nt = length(tv);

stab = dt/(dx^2);


kappa = 10;

xbar = -8.0;                      % Mitten av vågfunk. vid t=0
E = 0;
k0 =sqrt( E*2*m )/ hbar;                        % Medel vågtal/energi vid t=0, kopplad till r.mängd via p=hbar*k

sigmax = 0.4;                   % Standardavvikelse i x vid t=0
sigmak = 1/(2*sigmax);          % Standardavvikelse i k vid t=0

aa = 1/(4*sigmax^2);            % Konstant från uppgift 2.21 Griffiths
vg = hbar*k0/m;                 % Grupphastighet Verkar stämma!



% Normaliserat vågpaket som initialvillkor för numerisk lösare
Anorm = 1/(2*pi*sigmax^2)^(1/4);




Ux = 0.1*kappa*(xv./2).^4 - 2*kappa*(xv./2).^2;

H = @(x,k) (k.^2)/2 + 0.1*kappa*(x./2).^4 - 2*kappa*(x./2).^2 -100;

xtemp = -30:0.01:30;
ktemp = -40:0.1:40;
[xtest, ktest] = meshgrid(xtemp, ktemp);

Htest = H(xtest, ktest);


[egenfuncs, egenv, conv] = egenfunc(Nx+1, Ux, m, hbar, dx, 8);

conv

psi0 = Anorm*exp(-aa*(xv-xbar).^2) .* exp(1i*k0*(xv));

figure(2)
plot(xv, 50*egenfuncs(:,2))
hold on
plot(xv, 10*abs(psi0).^2, xv, Ux, '--')
yline(egenv(1,1))
yline(egenv(2,2))
yline(egenv(3,3))
hold off
ylim([-180, 100])
xlim([-12,12])



normtest = trapz(xv, abs(psi0).^2);


k0_v = fftshift( (1/sqrt(2*pi))*fft(psi0(1:Nx))*dx );
dk = 2*pi/(dx*Nx);
kv = ((-Nx/2:Nx/2-1)*dk).';
Nk = length(kv);

figure(3)
plot(kv, abs(k0_v).^2)

figure(4)
contour(xtemp,ktemp, Htest,-80:40:400)

T_o = sqrt(m/kappa)*2*pi;
%%
CN = CN_solve(psi0,Ux,Nx,Nt,dx,dt,m);
%%
max1 = max(abs(CN(1,:)));
max2 = max(abs(CN(end,:)));
%%
close all
nf = 100;
skip_frame = floor((Nt-1) / nf);
figure('Position',[100, 100, 1000, 900])
F = struct('cdata', cell(1,nf), 'colormap', cell(1,nf));

s1 = subplot(2,1,1);
s1.Subtitle.Interpreter = 'latex';
s1.Subtitle.FontSize = 12;
plot(xv, Ux, '-k')
hold on
PSI = plot(xv, 50*abs( CN(:,1) ).^2  );
xlabel('$x \ [\mathrm{eV}^{-1} \cdot \hbar c]$','fontsize',15,'Interpreter','latex')
legend('$U(x)$','$50|\Psi(t)|^2$','fontsize',15,'Interpreter','latex')
title('$U(x)=0.1\kappa(x/2)^4 - 2\kappa(x/2)^2, \ \kappa=10$ ','fontsize',18,'Interpreter','latex')
ylim([-100, 100])
xlim([-16,16])

s2 = subplot(2,1,2);
PHI = plot(kv(1:Nx), abs(fftshift( (1/sqrt(2*pi))*fft(CN(1:Nx, 1))*dx )));
xlabel('$k \ [\mathrm{eV} \cdot \hbar^{-1}c^{-1} ]$','fontsize',15,'Interpreter','latex')
ylabel('$|\Phi(k)|^2$','fontsize',15,'Interpreter','latex')
xlim([-50,50])
ylim([-0,1])

for j = 1:nf
    tj = sprintf('%.4f',skip_frame*j*dt);
    
    subtitle(s1, ['$t=$', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'],Interpreter="latex");
    
    PSI.YData = 50*abs( CN(:,skip_frame*j) ).^2  ;
    PHI.YData = abs(fftshift( (1/sqrt(2*pi))*fft(CN(1:Nx, j*skip_frame))*dx ));
    
    
    
%     hold off
    
    
    
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

figure('Position',[10, -10, 1600, 900])
F2 = struct('cdata', cell(1,nf), 'colormap', cell(1,nf));

zmax = max(max(max(Wv)));
zmin = min(min(min(Wv)));

[~,Cz] = contour( xtemp, ktemp, 100+Htest,-80:40:200,'-b'   );
Cz.ZLocation = 'zmax';
hold on
Wj = Wv(1:skip_l:end,1:skip_l:end,1).';
suf = surf(xv(1:skip_l:Nx), kv(1:skip_l:Nx), Wj(:,:,1));
suf.FaceColor="interp";
suf.EdgeColor = "none";
ylabel('$k \ [\mathrm{eV} \cdot \hbar^{-1}c^{-1} ]$','fontsize',15,'Interpreter','latex')
xlabel('$x \ [\mathrm{eV}^{-1} \cdot \hbar c]$','fontsize',15,'Interpreter','latex')

title('$\mathrm{Wigner-function}, \ U(x)=0.1\kappa(x/2)^4 - 2\kappa(x/2)^2, \ \kappa=10$', ...
    'fontsize',18,'Interpreter','latex')
xlim([-16,16])
ylim([-30,30])
%     zlim([zmin,20*zmax])
zlim([-100,400])
view(2)
colormap(flipud(cbrewer2('RdGy')));
clim([-5, 40])

c = colorbar;
ylabel(c,'$W(x,k)$','fontsize',15,'Interpreter','latex');
for j = 1:nf
    
    Wj = Wv(1:skip_l:end,1:skip_l:end,j).';
    suf.ZData = Wj;

    tj = sprintf('%.4f',skip_frame*j*dt);

    subtitle(['$t=$', num2str(tj), ' $[\hbar/ \mathrm{eV}]$'], Interpreter='latex', FontSize=12)
%         view(10, 45) 
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

psin = egenfuncs(1:Nx,nn);
En = egenv(nn,nn);

psim = egenfuncs(1:Nx,mm);
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
video=VideoWriter('doublewell_Wigner_2D_V2');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,F2);
close(video);

%%
video=VideoWriter('doublewell_psi_phi_V2','MPEG-4');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,F);
close(video);
