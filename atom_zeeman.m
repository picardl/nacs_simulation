
% atomic zeeman shifts for Na+Cs

clear;

c = constants();

Na = atom_basis('Na');
Cs = atom_basis('Cs');

Na = couple_angmom(Na,'s_Na','i_Na','f_Na');
Cs = couple_angmom(Cs,'s_Cs','i_Cs','f_Cs');

NaCs = join_basis(Na,Cs);

B1 = linspace(0,1000e-4,1e2);
[~,D] = eigenshuffle(NaCs.ops.H_Na_hyperfine+NaCs.ops.H_Cs_hyperfine...
    +(NaCs.ops.H0_Cs_zeeman+NaCs.ops.H0_Na_zeeman).*reshape(B1,1,1,[]));

figure(1); clf;
plot(B1*1e4,D*1e-9/c.h);
set(gca,'fontsize',14);
xlabel('B (Gauss)');
ylabel('E (GHz)');