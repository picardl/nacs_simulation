Na = DLine_basis('Na','D2');
c = constants();

Na = couple_angmom(Na,'j_NaD1','i_NaD1','f_NaD1');

[V,D] = eig(Na.ops.H_NaD2_hyperfine + Na.ops.H0_NaD2_zeeman.*B);
allEs = sort(squeeze(real(diag(D))));

if cpl
    [~,i1] = evec_ind({'f_NaD2','m_f_NaD2'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
else %If in uncoupled basis, use mJ and mI
    [~,i1] = evec_ind({'m_j_NaD2','m_i_NaD2'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
end
