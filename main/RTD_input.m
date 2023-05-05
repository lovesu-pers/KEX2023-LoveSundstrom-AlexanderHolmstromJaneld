%% Input file
% Creators: Love Sundström, Alexander Holmström Janeld
% Last modified: 04/05-2023

% Absolute Temprature
T = 300; 

% Number of energy points
NE = 80;

% RTD 1 dimensions [nm]
d = 1.5;                    % Barrier
w = 4.5;                    % Well 
LE = 0;                     % Emitter spacer
LC = 0;                     % Collector spacer
Ld = 10;                    % Extensions 
Ld_dop = 20/2;              % Doped region
L = LE + 2*d + w + LC;      % Total length

%% RTD 2 dimensions [nm] 
% d = 1.4;    
% w = 5.1;
% LE = 0;%1.5;
% LC = 0;%1.5;
% Ld = 10;
% Ld_dop = 60; 
% L = LE + 2*d + w + LC;


% Refrence energy set to emitter conduction band edge [eV]
E_0 = 1.42;

x = 0.26; % RTD 1
% x = 0.30;  % RTD 2 

% Conductionband discontinuty GaAs->AlGaAs [eV]
U0 = 0.6*1.247*x; 

% Start of first barrier, and limits of space domain
xb = 0;
x0 = -214.3;
xfin = (2*d+2+LC+Ld)+1+200;

% Relative effective mass [-]
mw = 0.067;                  %GaAs
md = (0.063+0.083*x);        % Al_xGaAs


% Emitter/collector doping [cm^-3]
n3d = 1e18;   % RTD1
% n3d = 2e17;     % RTD2


%% Scaleing for calculations in Hartree atomic units
m_e = 1;
hbar = 1;
q = 1;

q_SI = 1.602176634e-19;
a0_SI = 5.291772109e-11; 
m0_SI = 9.1093837e-31;
nm_a0_sf = 1e-9/a0_SI ;
Eh_SI = 4.35974472220712e-18;
Eh_eV = Eh_SI/q_SI;
V_Hartree_SI = 27.2113862456;
hbar_SI = 1.054571817e-34;
j_sf = q_SI*Eh_SI/(hbar_SI*a0_SI^2);
c_sf = 1/(a0_SI)^3;
t_sf = hbar_SI/Eh_SI;

d = d*nm_a0_sf;    
w = w*nm_a0_sf;
LE = LE*nm_a0_sf;
LC = LC*nm_a0_sf;
Ld = Ld*nm_a0_sf;
Ld_dop = Ld_dop *nm_a0_sf;

kB =  3.167e-6;
beta = kB*T;

%% Energy calcualtions

% Eq. 2d doping 
n2d = (n3d*a0_SI^3/(pi*1e-6)) * Ld_dop ;        % a0^-2;

% Fermi energies 
E_0 = E_0/Eh_eV;
E_FE0 = pi*hbar^2*n2d/mw + E_0;         % Davies p.33 
E_FE = beta*log(exp(E_FE0/beta)-1);


U0 = U0/Eh_eV; 

% Energy vector
EvE = linspace(E_0+beta*1e-5,E_0+5*beta,NE);

%% Scaling of space domain

x0 = x0*nm_a0_sf;
xfin = xfin*nm_a0_sf;

L = L*nm_a0_sf;