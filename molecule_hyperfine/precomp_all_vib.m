r = constants;

r.B = 8.8*1e-4; %B Field in Gauss. updated FR in model is at 864.1
r.initState = 2; %1 for FB mol, 2 for atoms
r.basis = 'aFC';
r.recompute = 2;
r.Jmax = 3;%Max J to inlcude in rotational basis
r.mtot = [3:6]; %Range of mF to include in hyperfine basis
r.pol_up = 'sigpm';
r.pol_frac = 0.9;
r.c3Sigma.Be = 1.14*1e9*r.h; %10.1103/PhysRevA.105.063322
r.c3Sigma.alpha1 =  0.228*1e9*r.h; %0.340*1e9*c.h;
r.c3Sigma.alpha2 =  0.143*1e9*r.h;
r.c3Sigma.Gamma = 150e6*r.h;
r.power = 100e-3;
r.waist = 13.4e-6/2;
r.f_vib = 288694e9;
vib_data = load('../data/cbB_210703_014725.mat');
E_vib = vib_data.out.E*r.hartree + r.Cs_D12_weighted;
ind = find(abs(E_vib-r.f_vib*r.h) == min(abs(E_vib-r.f_vib*r.h)));
r.E_vib = E_vib(ind);
r.N_tot = 0;

prev0 = precomp_hp(r);

r.f_vib = 306490e9;
prev12 = precomp_hp(r);

r.B = 860*1e-4;
r.f_vib = 310630e9;
prev15 = precomp_ops(r);

r.f_vib = 311991e9;
prev16 = precomp_ops(r);

r.f_vib = 313348e9;
prev17 = precomp_ops(r);

r.f_vib = 314698e9;
prev18 = precomp_ops(r);

r.f_vib = 320010e9;
prev22 = precomp_ops(r);

r.f_vib = 323873e9;
prev25 = precomp_ops(r);

r.f_vib = 325129e9;
prev26 = precomp_ops(r);

fn = '../data/precomp_all_v.mat';
save(fn,'prev0','prev12','prev15','prev16','prev17','prev18','prev22','prev25','prev26')