function basis_troubleshooting

const = constants();

B = 855e-4;

N_tot = 2;

%% uncoupled basis
Nmin = 0;
Nmax = N_tot;

basis.UC.qnums = build_basis({'s_1','s_2'},{1/2,1/2},[1,1]);
basis.UC.ops = build_operators(basis.UC.qnums);
basis_mol.qnums = build_basis({'eta','Lambda','N'},{1,0,Nmin:Nmax},[0,0,1],'b');
basis_mol.ops = build_operators(basis_mol.qnums);
basis.UC = join_basis(basis.UC,basis_mol);

% spin-coupled basis
[basis.SC,basis.UC,basis.change.UC_SC] = couple_angmom(basis.UC,'s_1','s_2','S');
basis = rmfield(basis,'UC'); % don't need spin-uncoupled basis anymore

%% projectors onto singlet and triplet subspaces
basis.SC.ops.Proj_singlet = 1*diag(basis.SC.qnums.S==0);
basis.SC.ops.Proj_triplet = 1*diag(basis.SC.qnums.S==1);

% out = sphten_compose(2,basis.SC.ops.S,T2);
% S = cat(3,basis.SC.ops.S_m,basis.SC.ops.S_z,basis.SC.ops.S_p);
% basis.SC.ops.S_zmol = sphdot(S,nz);
% Szmol = operator_matrix(@wigD_element,{basis.b.qnums,basis.a.qnums},{'S','m_S','Sigma'},'S',-1:1,0,'S','m_S','Sigma');

basis.SC.ops.F_z = basis.SC.ops.N_z + basis.SC.ops.S_z;

%% form hund's case a & b bases
[basis.b,basis.SC,basis.change.SC_b] = couple_angmom(basis.SC,'N','S','J');

% transformation to hund's case a
Jmax = max(basis.b.qnums.J);
Jmin = 0;
basis.a.qnums = build_basis({'eta','s_1','s_2','J','S','Lambda'},...
    {1,const.s_1,const.s_2,Jmin:Jmax,0:1,0},[0 0 0 0 0 0],'a');

% del = (basis.a.qnums.J>abs(basis.a.qnums.Omega));
% basis.a.qnums(del,:) = [];

basis.a.ops = build_operators(basis.a.qnums);
basis.change.a_b = operator_matrix(@case_b2a_element,{basis.a.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});

% nz = operator_matrix(@wigD_element,{basis.b.qnums,basis.a.qnums},{'N','Lambda','m_N','S','m_S','Sigma'},'S','m_S','Sigma','N','Lambda','m_N');
nz = operator_matrix(@wigD_element,basis.SC.qnums,{'N','Lambda','m_N'},1,-1:1,0,'N','Lambda','m_N');
S = cat(3,basis.SC.ops.S_m,basis.SC.ops.S_z,basis.SC.ops.S_p);
basis.SC.ops.S_zmol = sphdot(S,nz);

end

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
c.apol_2 = 1.91395507e-38; % C*m^2/V

%% NaCs specific
c.i_2 = 0; %7/2;
c.s_2 = 1/2;
c.zeta_hf_2 = c.h*2.2981579425e9; % Cs hyperfine constant, Joules
c.gs_2 = 2.00254032; % g_S, spin g-factor in Cs ground state, using steck definition
c.gi_2 = -0.00039885395; % g_I, nuclear spin g-factor in Cs ground state, using steck definition

c.i_1 = 0; %3/2;
c.s_1 = 1/2;
c.zeta_hf_1 = c.h*885.81306440e6;
c.gs_1 = 2.0022960;
c.gi_1 = -0.00080461080;

c.m_nacs = (133+23)*c.amu;
c.mu_nacs = (133*23/(133+23))*c.amu;

c.c3Sigma.Lambda = 0;
c.c3Sigma.Omega = 1;
c.c3Sigma.S = 1;
% c.c3Sigma.Be = 0.925*1e9*c.h;
c.c3Sigma.Be = 0.962*1e9*c.h;
% c.c3Sigma.Be = 0.995*1e9*c.h;
c.c3Sigma.gS = 2.0023;
c.c3Sigma.alpha1 = 0.340*1e9*c.h;
% c.c3Sigma.alpha1 = -0.200*1e9*c.h;
c.c3Sigma.alpha2 = 0.095*1e9*c.h;
% c.c3Sigma.alpha2 = 0.250*1e9*c.h;
c.c3Sigma.wef = c.h*5e6;
c.c3Sigma.singlet_fraction = 0.322;
c.c3Sigma.Gamma = c.h*120e6;

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

c.Cs_D2 = c.h*351725718500813;
c.Cs_D1 = c.h*335116048808294;
c.Cs_D12_weighted = (4*c.Cs_D2 + 2*c.Cs_D1)/6;

end