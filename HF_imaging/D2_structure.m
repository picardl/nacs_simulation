basis = D2_basis('CsD2');
c = constants();

Cs = couple_angmom(basis,'j_CsD2','i_CsD2','f_CsD2');
[V,D] = eig(Cs.ops.H_CsD2_hyperfine);
out = sort(squeeze(real(diag(D))));
zeropt = min(out);
out = (out - zeropt)/c.h/1e6;
mfs = [-2:2,-3:3,-4:4,-5:5];
scatter(mfs,out);
xlabel('mF')
ylabel('Energy [MHz]')
title('Cs D2 HF structure, zero field')

[~,i1] = evec_ind({'f_CsD2','m_f_CsD2'},[5,5],Cs,V);
disp((real(D(i1,i1)) - zeropt)/c.h/1e6);