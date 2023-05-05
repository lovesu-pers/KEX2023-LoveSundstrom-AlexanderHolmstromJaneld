clear

RTD_input;
md = md;
w = 1*w;
d = 1.5*d;
Nx = 1500;
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


U = 1*U0;
[Ux,~,mx ,~,~ ,~,NR,~,~,~] = ... 
                        simpUx_mxV5(xv,(U),Vb,xb,w,d,LE,LC,q,md,mw,Ec);

xbar = -1000;

Ev = linspace(0.001*U,2*U,80); %E_FE-Ec;

%



% kspace = fftshift(1/(sqrt(2*pi)) * fft(psi0)*dx);
% dk = 2*pi/(dx*Nx);
% kv = (-Nx/2 : Nx/2-1)*dk;
% 
%  (Nx/2-1)*dk

int_index = (Nx-NR)
% 
% figure(1)
% yyaxis left
% plot(Ux)
% yyaxis right
% plot(abs(psi0) )
% ylim([0,0.50])
% xlim([550 1000])
% 
% 
% figure(2)
% plot(kv,abs(kspace).^2)
%%

for j = 1:length(Ev)
    E0 = Ev(j);
    k0 =  sqrt(2*mw*(E0) )/(hbar);

    vg1 =  hbar*sqrt(2*mw*(E0) )/(hbar*mw);
    vg2 =  hbar*sqrt(2*mw*(E0-Ec) )/(hbar*mw);
    
    
    
    tmax = ceil((2000-xbar)/(vg1));
    
    % dt = tmax/(Nt);
    tv = linspace(0,tmax,Nt);
    dt = tv(2)-tv(1);
    
    % vg2*tmax
    
    sigmax = 4*d;
    aa = 1/(4*sigmax^2); 
    
    Anorm = 1/(2*pi*sigmax^2)^(1/4);
    psi0 = Anorm*exp(-aa*(xv-xbar).^2) .* exp(1i*k0*(xv));
    normtest = trapz(xv, abs(psi0).^2);

PSI_CN = CN_solve_mx(psi0, Ux, mx, Nx, Nt, Nframe, dx, dt, hbar,0);
Tv(j) = trapz(xv(int_index:end), abs(PSI_CN(int_index:end,end)).^2);
j
end

%%
figure(1)
plot(Ev,Tv)