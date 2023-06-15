%% Model Setup
%Voltage Noise
noise_pow = 1e-07;
noise_pow_exp = 1e-13;%1e-14;
% Cathode Observer Gain (<1/4)
lambda = -20;
% Voltage Inversion Gain
gamma = 1e8;
% Anode Observer Gain
k_anode = 1e-2;
% exp inv gain
% gamma_exp = 1e3;
gamma_exp = 1e22;
%% Parameters
% Stoichiometric Windows 
CC.x100_percent_n = 0.8292;
CC.x0_percent_n = 3.5506e-04; 
CC.y0_percent_p = 0.8941;
CC.y100_percent_p = 0.0335; 
Cap_n =  6.0272;
Cap_p =  5.8050;
Capacity = 4.9872;

% Constants
CC.F = 96485;    % Faraday Constant (C/mol)
CC.R = 8.3145;   % Gas constant
CC.T = 273.15+25;% Temperature
CC.Tref = 273.15;

% Electrode Parameters
CC.Ds_n = 5e-15; % Ds_n (m^2/s)         Graphite 0.4
CC.Ds_p = 8e-15; % Ds_p_1 (m^2/s)      NMC 0.025
nr = 50;
CC.nr = nr; % Number of elements in the solid particle
CC.nn = 50;  %number of elements in the negative electrode
CC.ns = 25;  %number of elements in the separator
CC.np = 50; %number of elements in the positive electrode
CC.Rn_0 = 2.5e-06;  % Radius of Particle       Graphite
CC.Rp_0 = 3.5e-06; % Radius of Particle       NMC
CC.del_n = 6.2e-05; % Electrode Thickness         Graphite
CC.del_p = 6.7e-05; % Electrode Thickness        NMC
CC.del_s = 1.2e-05; % Separator Thickness

CC.drn_0 = CC.Rn_0/CC.nr;
CC.drp_0 = CC.Rp_0/CC.nr;   
CC.kt_n = 1*0.0250;   % kt (omega m^2) film resistence  Graphite
CC.kt_p = 0*300e-4;  % kt (omega m^2) film resistence  NMC
CC.k0_n = 1.061e-6/CC.F; % k0_n (mol/[m^2s(mol/m^3)]) reaction rate constant  Graphite
CC.k0_p = 4.824e-6/CC.F; % k0_n (mol/[m^2s(mol/m^3)]) reaction rate constant  NMC

CC.cs_max_n = 28746; % cs_max_n (mol/m^3)      Graphite
CC.cs_max_p = 35380; % cs_max_p_1 (mol/m^3)    NMC

% Electrolyte Parameters
CC.ce_init = 1000;     %Initial Electrolyte Concentration (mol/m^3)
CC.tf = 0.62;          % Thermodynamic factor (averaged) *(1-t_plus)*(1+d ln(f)/ d ln(c_e)) with c_e = 1 mol/l
CC.t_plus = 0.38;      % transference number
CC.kappa = 1.3; %1.4       % (S/m) electrolyte ionic conductivity (averaged)
CC.D_e = 5.35e-10;      % (m^2/s) Electrolyte Diffusion Coefficient
CC.brugg = 1.5;        % Bruggman exponent

CC.epsilon_e_n = 0.3;   % Volume fraction in electrolyte (porosity) for neg. electrode
CC.epsilon_e_s = 0.4;   % Volume fraction in electrolyte (porosity) for separator
CC.epsilon_e_p = 0.3;   % Volume fraction in electrolyte (porosity) for pos. electrode

CC.A = 0.205; %Battery Total Area (Area x No of Layers)
% Individual Electrode Capacity Calculation
% Cap_n = CC.F*CC.es_n*CC.A*CC.del_n*CC.cs_max_n/3600;
% Cap_p = CC.F*CC.es_p*CC.A*CC.del_p*CC.cs_max_p/3600;
es_n = Cap_n*3600/(CC.F*CC.A*CC.del_n*CC.cs_max_n);
es_p = Cap_p*3600/(CC.F*CC.A*CC.del_p*CC.cs_max_p);
CC.es_n = es_n;
CC.es_p = es_p;
CC.as_n_0 = 3*CC.es_n/CC.Rn_0; % Active surface area per electrode unit volume   Graphite
CC.as_p_0 = 3*CC.es_p/CC.Rp_0; % Active surface area per electrode unit volume   NMC

% Mechanical Properties
CC.E_nP = 15e9; % Young's Modulus
CC.v_nP = 0.3; % Poisson's Ratio
CC.Omega_nP = 4.08154e-6; % Partial Molar Volume
CC.k_exp = 96; % Tuning parameter for constrained expansion


% Thermal properties
CC.h = 0.9462; %The lumped convective coefficient W/m^2/K
CC.th_mass = 1864.54;% Thermal Mass J/m^3/K
CC.alp_t =  2e-08; % Expansion Thermal Coefficient
%% Arrhenius 
CC.E_Ds_p = 18550;
CC.E_r_p = 39570;
CC.E_Ds_n = 42770;
CC.E_r_n = 37480;
CC.E_De = 37040;
CC.E_ke = 34700;

%% SEI Plating parameters
CC.U_sei = 0;
CC.i0pl = 0.001;
CC.del_sei_init = 5e-9;
CC.D_ec = 2e-18;
CC.c_ec0 = 4.541;
CC.M_sei = 0.162;
CC.rho_sei = 1690;
CC.k0_sei = 1e-12;
CC.k_sei = 5e-6;
CC.M_li = 6.94e-3;
CC.rho_li = 534;


%% Particle Diffusion Matrices for Plant (Method of Lines)

Asys=zeros(nr+1,nr+1);
Asys(1,1)=-6;
Asys(1,2)=6;
for i=1:nr-1
    Asys(i+1,i)=(i-1)/i;
    Asys(i+1,i+1)=-2;
    Asys(i+1,i+2)=(i+1)/i;
end
Asys(nr+1,nr)=2;
Asys(nr+1,nr+1)=-2;
Bsys=zeros(nr+1,1);
Bsys(nr+1,1)=-2*(nr+1)/nr;

CC.Asys = Asys;
CC.Bsys = Bsys;


%% Electrolyte Equations FC
epsilon_e_s = CC.epsilon_e_s;
epsilon_e_p = CC.epsilon_e_p;
epsilon_e_n = CC.epsilon_e_n;
CC.dx_n = CC.del_n/CC.nn;
CC.dx_s = CC.del_s/CC.ns;
CC.dx_p = CC.del_p/CC.np;

dx_n = CC.dx_n;
dx_s = CC.dx_s;
dx_p = CC.dx_p;
np = CC.np;
ns = CC.ns;
nn = CC.nn;

CC.del_t = CC.del_n+CC.del_s+CC.del_p;
x_neg = 0:dx_n:CC.del_n;
x_sep = CC.del_n:dx_s:CC.del_n+CC.del_s;
x_pos = CC.del_n+CC.del_s:dx_p:CC.del_n+CC.del_s+CC.del_p;
x_bat = [x_neg(1:end-1) x_sep(1:end-1) x_pos]/CC.del_t;

De_s=CC.D_e*epsilon_e_s^CC.brugg;
De_p=CC.D_e*epsilon_e_p^CC.brugg;
De_n=CC.D_e*epsilon_e_n^CC.brugg;


% Positive | Seperator | Negative
A=diag([De_p/epsilon_e_p/dx_p/dx_p*(-2)*ones(1,np+1),De_s/epsilon_e_s/dx_s/dx_s*(-2)*ones(1,ns),De_n/epsilon_e_n/dx_n/dx_n*(-2)*ones(1,nn)])+...
    diag([De_p/epsilon_e_p/dx_p/dx_p*ones(1,np),De_s/epsilon_e_s/dx_s/dx_s*ones(1,ns),De_n/epsilon_e_n/dx_n/dx_n*ones(1,nn)],1)+...
    diag([De_p/epsilon_e_p/dx_p/dx_p*ones(1,np),De_s/epsilon_e_s/dx_s/dx_s*ones(1,ns),De_n/epsilon_e_n/dx_n/dx_n*ones(1,nn)],-1);

A(1,2)=2*De_p/epsilon_e_p/dx_p/dx_p;

A(np+1,np+1-1)=(2*De_p*dx_s )/(dx_p*dx_s*(dx_p*epsilon_e_p + dx_s*epsilon_e_s));% Cep[N-1]
A(np+1,np+1)=(-2*De_p*dx_s-2*De_s*dx_p    )/(dx_p*dx_s*(dx_p*epsilon_e_p + dx_s*epsilon_e_s));% Cep[N}= Ces[0]
A(np+1,np+1+1)=(2*De_s*dx_p  )/(dx_p*dx_s*(dx_p*epsilon_e_p + dx_s*epsilon_e_s)); % Ces[1]

A(np+ns+1,np+ns+1-1)=(2*De_s)/(dx_s*(dx_n*epsilon_e_n + dx_s*epsilon_e_s));% Cep[N-1]
A(np+ns+1,np+ns+1)=-(2*(De_n*dx_s + De_s*dx_n))/(dx_n*dx_s*(dx_n*epsilon_e_n + dx_s*epsilon_e_s));% Cep[N}= Ces[0]
A(np+ns+1,np+ns+1+1)=(2*De_n)/(epsilon_e_n*dx_n^2 + dx_s*epsilon_e_s*dx_n); % Ces[1]

A(np+ns+nn+1,np+ns+nn+1)=-2*De_n/epsilon_e_n/dx_n/dx_n;
A(np+ns+nn+1,np+ns+nn)=2*De_n/epsilon_e_n/dx_n/dx_n;

B=zeros(np+ns+nn+1,1);

II=1:(np+1);
kkp=(1-CC.t_plus)/epsilon_e_p/CC.F/CC.del_p;
B(II)=-kkp;% % const J
B(np+1)= -(dx_p^2*dx_s*epsilon_e_p*kkp)/(dx_p*dx_s*(dx_p*epsilon_e_p + dx_s*epsilon_e_s));

IL = np+ns+1:np+ns+nn+1;
kkn=(1-CC.t_plus)/epsilon_e_n/CC.F/CC.del_n;
B(IL)=kkn;% % const J
B(np+ns+1)= (dx_n^2*dx_s*epsilon_e_n*kkn)/(dx_n*dx_s*(dx_n*epsilon_e_n + dx_s*epsilon_e_s));

%C=eye(np+ns+nn+1)
C=[ones(1,np+1)*dx_p,ones(1,ns)*dx_s,ones(1,nn)*dx_n];
C(1)=dx_p/2;
C(np+1)=dx_p/2+dx_s/2;
C(np+ns+1)=dx_n/2+dx_s/2;
C(np+ns+nn+1)=dx_n/2;

CC.Ael = A;
CC.Bel = B;
CC.Cel = C;

CC.n_nodes = nn+ns+np+1;
%% Initialize Plant and Observer Concentrations
CC.SOC_0 = SOC_0;
% Initial Particle Concentrations (Plant)
cs_init_n =   CC.cs_max_n*(CC.SOC_0*(CC.x100_percent_n-CC.x0_percent_n)+CC.x0_percent_n);
% cs_init_n = 48.8682;
cs_init_p =   CC.cs_max_p*(CC.SOC_0*(CC.y100_percent_p-CC.y0_percent_p)+CC.y0_percent_p);
% cs_init_p = 31513.0;
% Total Lithium
nLi = CC.es_p*CC.del_p*cs_init_p+CC.es_n*CC.del_n*cs_init_n;




