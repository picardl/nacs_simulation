
clear;

c = constants();

Na = atom_basis('Na');
Cs = atom_basis('Cs');

Na = couple_angmom(Na,'s_Na','i_Na','f_Na');
Cs = couple_angmom(Cs,'s_Cs','i_Cs','f_Cs');

NaCs = join_basis(Na,Cs);

B1 = linspace(0,865e-4,1e3);
[~,D] = eigenshuffle(NaCs.ops.H_Na_hyperfine+NaCs.ops.H_Cs_hyperfine...
    +(NaCs.ops.H0_Cs_zeeman+NaCs.ops.H0_Na_zeeman).*reshape(B1,1,1,[]));
% [V_Cs,D_Cs] = eigenshuffle(Cs.ops.H_Cs_hyperfine+Cs.ops.H0_Cs_zeeman.*reshape(B1,1,1,[]));

% B2 = linspace(0,865e-4,1e3);
% [V_Na,D_Na] = eigenshuffle(Na.ops.H_Na_hyperfine+Na.ops.H0_Na_zeeman.*reshape(B2,1,1,[]));

figure(1); clf;
% subplot(2,1,1);
plot(B1,D*1e-9/c.h,'-b');
% ylim([-25 25])

% subplot(2,1,2);
% plot(B2,D_Na*1e-9/c.h,'-b');
% ylim([-5 5])