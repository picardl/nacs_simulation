function c = constants()

%% fundamental constants
c.alpha = 0.00729735256;
c.c = 299792458;
c.hbar = 1.054571817646156e-34;
c.me =9.1093837015e-31; % kg
c.h = 2*pi*c.hbar;
c.amu = 1.66053906660e-27;
c.eta0 = 376.730313; % ohms
c.eps0 = 8.85418782e-12;
c.uB = 9.2740100783e-24; % J/T
c.uN = 5.050783699e-27; % J/T
c.e = 1.6021766340e-19; % C
c.mu0 = 4*pi*1e-7;
c.kb = 1.380649e-23; %J/K

%% derived constants
c.abohr = 5.29177210903e-11;
c.hartree = c.hbar^2/(c.me*c.abohr^2);
c.B_au = c.hartree/(2*c.uB);
c.wavenum2hartree = 4.55634011e-6;
c.apol_Cs = 1.91395507e-38; % C*m^2/V

%% NaCs specific
c.i_Cs =  7/2;
c.s_Cs = 1/2;
c.zeta_hf_Cs = c.h*2.2981579425e9; % Cs hyperfine constant, Joules
c.gs_Cs = 2.00254032; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_Cs = -0.00039885395; % g_I, nuclear spin g-factor in Cs ground state, using steck definition

c.i_CsD2 =  7/2;
c.s_CsD2 = 1/2;
c.j_CsD2 = 3/2;
c.zeta_hf_CsD2 = c.h*50.28827e6; % Cs hyperfine constant, Joules
c.B_hf_CsD2 = -1*c.h*0.4934e6; % Cs hyperfine constant, Joules
c.C_hf_CsD2 = c.h*0.560e3; % Cs hyperfine constant, Joules
c.gj_CsD2 = 1.33400; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_CsD2 = -0.00039885395; % g_I, nuclear spin g-factor in Cs ground state, using steck definition
c.gs_CsD2 = 2.00254032;

c.i_CsD1 =  7/2;
c.s_CsD1 = 1/2;
c.j_CsD1 = 1/2;
c.zeta_hf_CsD1 = c.h*291.9201e6; % Cs hyperfine constant, Joules
c.B_hf_CsD1 = 0; % Cs hyperfine constant, Joules
c.C_hf_CsD1 = 0; % Cs hyperfine constant, Joules
c.gj_CsD1 = 0.665900; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_CsD1 = -0.00039885395; % g_I, nuclear spin g-factor in Cs ground state, using steck definition
c.gs_CsD1 = 2.00254032;

c.i_Na = 3/2;
c.s_Na = 1/2;
c.zeta_hf_Na = c.h*885.81306440e6;
c.gs_Na = 2.0022960;
c.gi_Na = -0.00080461080;


c.i_NaD1 =  3/2;
c.s_NaD1 = 1/2;
c.j_NaD1 = 1/2;
c.zeta_hf_NaD1 = c.h*94.44e6; % Na hyperfine constant, Joules
c.B_hf_NaD1 = 0; % Cs hyperfine constant, Joules
c.C_hf_NaD1 = 0; % Cs hyperfine constant, Joules
c.gj_NaD1 = 0.66581; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_NaD1 = -0.00080461080; % g_I, nuclear spin g-factor in Cs ground state, using steck definition
c.gs_NaD1 = 2.0023193043622;

%energy differences in Hz
c.dE_NaD2 = (-16973.36619+[0
    25739.999
    29172.837
    29172.887
    33200.673
    34548.764
    36372.618])*29.9792458e9;
c.dE_NaD1 = (-16956.17025+[0
    25739.999
    29172.887
    33200.673
    34548.764
    36372.618])*29.9792458e9;
c.dE_Na = ([16956.17025
    16973.36619 
    25739.999
    29172.887
    33200.673
    34548.764
    36372.618])*29.9792458e9;

c.dEj_NaD2 = [1/2,1/2,5/2,3/2,1/2,3/2,5/2,1/2];
c.dEj_NaD1 = [1/2,1/2,3/2,1/2,3/2,5/2,1/2];
c.dEj_Na = [1/2,3/2,1/2,3/2,1/2,3/2,5/2,1/2];

function atomPol(a,J,K,dE,dEj)
   (-1).^(K + J + 1).*sqrt(2*K+1) 
end

c.i_NaD2 =  3/2;
c.s_NaD2 = 1/2;
c.j_NaD2 = 1/2;
c.zeta_hf_NaD2 = c.h*18.534e6; % Na hyperfine constant, Joules
c.B_hf_NaD2 = c.h*2.724e6; % Cs hyperfine constant, Joules
c.C_hf_NaD2 = 0; % Cs hyperfine constant, Joules
c.gj_NaD2 = 0.66581; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_NaD2 = -0.00080461080; % g_I, nuclear spin g-factor in Cs ground state, using steck definition
c.gs_NaD2 = 2.0023193043622;

c.m_nacs = (132.905451933+22.9897692820)*c.amu;
c.mu_nacs = (132.905451933*22.9897692820/(132.905451933+22.9897692820))*c.amu;

c.c3Sigma.Lambda = 0;
c.c3Sigma.Omega = 1;
c.c3Sigma.S = 1;
% c.c3Sigma.Be = 0.925*1e9*c.h;
% c.c3Sigma.Be = 0.962*1e9*c.h;
% c.c3Sigma.Be = 0.995*1e9*c.h;
c.c3Sigma.Be = 1.139*1e9*c.h; %10.1103/PhysRevA.105.063322
c.c3Sigma.alpha1 = 0.4429*1e9*c.h; %0.340*1e9*c.h;
c.c3Sigma.alpha2 = -0.09568*1e9*c.h;
c.c3Sigma.gamma = 2*1.139*1e9*c.h;
c.c3Sigma.lambda = 0*40*1e9*c.h;
c.c3Sigma.gS = 2.0023;
% c.c3Sigma.alpha1 = -0.200*1e9*c.h;
% c.c3Sigma.alpha2 = 0.095*1e9*c.h;
c.c3Sigma.lI = 0.085549*1e9*c.h;
% c.c3Sigma.alpha2 = 0.250*1e9*c.h;
c.c3Sigma.wef = c.h*15e6;
c.c3Sigma.singlet_fraction = 0.075;
c.c3Sigma.Gamma = c.h*15e6;

c.X1Sigma.Bv = 1.7356e9*c.h; % rot const
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
c.X1Sigma.a_perp = 1872.1153;
c.X1Sigma.a_par = 467.0375; 

c.Cs_D2 = c.h*351725718500813;
c.Cs_D1 = c.h*335116048808294;
c.Cs_D12_weighted = (4*c.Cs_D2 + 2*c.Cs_D1)/6;

end