%Calculate and plot Zeeman shifts of Na ground, D1 or D2 lines
c = constants();

BFields = [0:0.1:8]*1e-4;
elevs = zeros(8,0);
Jm_Na = [2,2];
cpl = 0;

for i = BFields
%     [allE,V,eT] = Na_HF_B(i);%,cpl,Jm_Na);
    allE = Na_HF_B(i);%,cpl,Jm_Na);
    elevs(:,end+1)  = allE;
end

plot(BFields*1e4,(elevs-min(elevs,[],1))/c.h/1e6);
xlabel('B Field [G]')
ylabel('Energy [MHz]')
title('Na D2 Zeeman')
disp(['Energy of [',num2str(Jm_Na(1)),',',num2str(Jm_Na(2)),'] state = ',num2str(eT/c.h/1e6),' MHz']);
function out = Na_HF_B(B)

if nargin<1
    fm_Na = [4 4];
end

Na = atom_basis('Na');

Na = couple_angmom(Na,'s_Na','i_Na','f_Na');

[V,D] = eig(Na.ops.H_Na_hyperfine + Na.ops.H0_Na_zeeman.*B);
out = sort(squeeze(real(diag(D))));

end

function [allEs,V,targetE] = NaD2_HF_B(B,cpl,fm_Na)

if nargin<2
    cpl = 1;
end
if nargin < 3 && cpl == 1;
    fm_Na = [5,5];
end

Na = DLine_basis('Na','D2');

if cpl %Boolean, whether or not to couple J and I
    Na = couple_angmom(Na,'j_NaD2','i_NaD2','f_NaD2');
end

[V,D] = eig(Na.ops.H_NaD2_hyperfine + Na.ops.H0_NaD2_zeeman.*B);
allEs = sort(squeeze(real(diag(D))));

if cpl
    [~,i1] = evec_ind({'f_NaD2','m_f_NaD2'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
else %If in uncoupled basis, use mJ and mI
    [~,i1] = evec_ind({'m_j_NaD2','m_i_NaD2'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
end

end

function [allEs,V,targetE] = NaD1_HF_B(B,cpl,fm_Na)

if nargin<2
    cpl = 1;
end
if nargin < 3 && cpl == 1
    fm_Na = [4,4];
end

Na = DLine_basis('Na','D1');

if cpl %Boolean, whether or not to couple J and I
    Na = couple_angmom(Na,'j_NaD1','i_NaD1','f_NaD1');
end

[V,D] = eig(Na.ops.H_NaD2_hyperfine + Na.ops.H0_NaD2_zeeman.*B);
allEs = sort(squeeze(real(diag(D))));

if cpl
    [~,i1] = evec_ind({'f_NaD1','m_f_NaD1'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
else %If in uncoupled basis, use mJ and mI
    [~,i1] = evec_ind({'m_j_NaD2','m_i_NaD2'},fm_Na,Na,V);
    targetE = real(D(i1,i1));
end

end

