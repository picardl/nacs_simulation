%Calculate and plot Zeeman shifts of Cs ground, D1 or D2 lines
c = constants();

BFields = [0:50:850]*1e-4;
elevs = zeros(16,0);
Jm_Cs = [-3/2,7/2];
cpl = 0;

for i = BFields
    [allE,V,eT,targeteV] = Cs_HF_B(i,[4,4]);
    elevs(:,end+1)  = allE;
end

plot(BFields*1e4,elevs/c.h/1e6);
xlabel('B Field [G]')
ylabel('Energy [MHz]')
title('Cs D2 Zeeman')
disp(['Energy of [',num2str(Jm_Cs(1)),',',num2str(Jm_Cs(2)),'] state = ',num2str(eT/c.h/1e6),' MHz']);

function [allEs,V,targetE,targeteV] = Cs_HF_B(B, fm_Cs)

if nargin<1
    fm_Cs = [4 4];
end

Cs = atom_basis('Cs');

Cs = couple_angmom(Cs,'s_Cs','i_Cs','f_Cs');

[V,D] = eig(Cs.ops.H_Cs_hyperfine + Cs.ops.H0_Cs_zeeman.*B);
[allEs,sortInd] = sort(squeeze(real(diag(D))));
[~,i1] = evec_ind({'f_Cs','m_f_Cs'},fm_Cs,Cs,V);
targetE = real(D(i1,i1));
targeteV = V(:,i1);  
end

function [allEs,V,targetE,targeteV] = CsD2_HF_B(B,cpl,fm_Cs)

if nargin<2
    cpl = 1;
end
if nargin < 3 && cpl == 1
    fm_Cs = [5,5];
end

Cs = DLine_basis('Cs','D2');

if cpl %Boolean, whether or not to couple J and I
    Cs = couple_angmom(Cs,'j_CsD2','i_CsD2','f_CsD2');
end

[V,D] = eig(Cs.ops.H_CsD2_hyperfine + Cs.ops.H0_CsD2_zeeman.*B);
allEs = sort(squeeze(real(diag(D))));

if cpl
    [~,i1] = evec_ind({'f_CsD2','m_f_CsD2'},fm_Cs,Cs,V);
    targetE = real(D(i1,i1));
else %If in uncoupled basis, use mJ and mI
    [~,i1] = evec_ind({'m_j_CsD2','m_i_CsD2'},fm_Cs,Cs,V);
    targetE = real(D(i1,i1));
end
targeteV = V(:,i1);  
end

function [allEs,V,targetE,targeteV] = CsD1_HF_B(B,cpl,fm_Cs)

if nargin<2
    cpl = 1;
end
if nargin < 3 && cpl == 1
    fm_Cs = [5,5];
end

Cs = DLine_basis('Cs','D1');

if cpl %Boolean, whether or not to couple J and I
    Cs = couple_angmom(Cs,'j_CsD1','i_CsD1','f_CsD1');
end

[V,D] = eig(Cs.ops.H_CsD1_hyperfine + Cs.ops.H0_CsD1_zeeman.*B);
allEs = sort(squeeze(real(diag(D))));

if cpl
    [~,i1] = evec_ind({'f_CsD1','m_f_CsD1'},fm_Cs,Cs,V);
    targetE = real(D(i1,i1));
else %If in uncoupled basis, use mJ and mI
    [~,i1] = evec_ind({'m_j_CsD1','m_i_CsD1'},fm_Cs,Cs,V);
    targetE = real(D(i1,i1));
end
targeteV = V(:,i1);  
end