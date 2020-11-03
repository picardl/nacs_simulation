function c = constants()

%% fundamental constants
c.alpha = 0.00729735256;
c.c = 299792458;
c.hbar = 1.0545718e-34;
c.me = 9.10938356e-31; % kg
c.h = 2*pi*c.hbar;
c.amu = 1.66054e-27;
c.eta0 = 376.730313; % ohms
c.eps0 = 8.85418782e-12;
c.uB = 9.274009994e-24; % J/T
c.uN = 5.050783699e-27; % J/T
c.e = 1.60217662e-19; % C
c.mu0 = 4*pi*1e-7;

%% derived constants
c.abohr = 5.29177210903e-11;
c.hartree = c.hbar^2/(c.me*c.abohr^2);
c.B_au = c.hartree/(2*c.uB);
c.wavenum2hartree = 4.55634011e-6;
c.apol_Cs = 1.91395507e-38; % C*m^2/V

%% NaCs specific
c.i_Cs = 7/2;
c.s_Cs = 1/2;
c.zeta_hf_Cs = c.h*2.2981579425e9; % Cs hyperfine constant, Joules
c.gs_Cs = 2.0023193043737; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_Cs = -0.00039885395; % g_I, nuclear spin g-factor in Cs ground state, using steck definition

c.i_Na = 3/2;
c.s_Na = 1/2;
c.zeta_hf_Na = c.h*885.81306440e6;
c.gs_Na = 2.0023193043622;
c.gi_Na = -0.00080461080;

c.m_nacs = (133+23)*c.amu;
c.mu_nacs = (133*23/(133+23))*c.amu;

c.c3Sigma.Lambda = 0;
c.c3Sigma.Omega = 1;
c.c3Sigma.S = 1;
c.c3Sigma.Be = 0.925*1e9*c.h;
c.c3Sigma.gamma = 2*c.c3Sigma.Be;
c.c3Sigma.gS = 2.0023;
c.c3Sigma.alpha1 = 0.340*1e9*c.h;
c.c3Sigma.alpha2 = 0.095*1e9*c.h;
c.c3Sigma.wef = c.h*5e6;
c.c3Sigma.singlet_fraction = 0.5;
c.c3Sigma.Gamma = c.h*90e6;

c.X1Sigma.Bv = 1.7396e9*c.h; % rot const
c.X1Sigma.eQq1 = -0.097e6*c.h; % Na elec quadrupole const
c.X1Sigma.eQq2 = 0.150e6*c.h; % Cs elec quadrupole const
c.X1Sigma.c1 = 14.2*c.h;
c.X1Sigma.c2 = 854.5*c.h;
c.X1Sigma.c3 = 105.6*c.h;
c.X1Sigma.c4 = 3941.8*c.h;
c.X1Sigma.gr = 1e-3; % rotational g-factor
c.X1Sigma.g1 = 1.478; % Na nuc g-factor
c.X1Sigma.g2 = 0.738; % Cs nuc g-factor
c.X1Sigma.sigma1 = 639.2e-6; % Na nuc spin shielding
c.X1Sigma.sigma2 = 6278.7e-6; % Cs nuc spin shielding

end