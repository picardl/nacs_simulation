
function out = atom_zeeman(fm_Cs,fm_Na,B)

if nargin<3
    B = 10e-4;
end
if nargin<2
    fm_Na = [1 1];
end
if nargin<1
    fm_Cs = [3 3];
end

Na = atom_basis('Na');
Cs = atom_basis('Cs');

Na = couple_angmom(Na,'s_Na','i_Na','f_Na');
Cs = couple_angmom(Cs,'s_Cs','i_Cs','f_Cs');

NaCs = join_basis(Na,Cs);

[V,D] = eig(NaCs.ops.H_Na_hyperfine+NaCs.ops.H_Cs_hyperfine...
    +(NaCs.ops.H0_Cs_zeeman+NaCs.ops.H0_Na_zeeman).*B);

[~,i1] = evec_ind({'f_Cs','m_f_Cs','f_Na','m_f_Na'},[fm_Cs,fm_Na],NaCs,V);
out = real(D(i1,i1));

end
