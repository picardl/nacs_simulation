%%Na
B = 860*1e-4;
c = constants();

%Generate bases and Hamiltonians, both uncoupled (mJ,mI) and coupled (F,mF)
%bases
NaGnd = atom_basis('Na');
NaD2 = DLine_basis('Na','D2');
NaGnd_cpl = couple_angmom(NaGnd,'s_Na','i_Na','f_Na');
NaD2_cpl = couple_angmom(NaD2,'j_NaD2','i_NaD2','f_NaD2');

%% Define reference point of 0 magnetic field transition
ref_g = [2,2]; %[F,mF] in coupled basis
ref_e = [3,3]; %[F,mF] in coupled basis

[V,D] = eig(NaGnd_cpl.ops.H_Na_hyperfine);
[~,i1] = evec_ind({'f_Na','m_f_Na'},ref_g,NaGnd_cpl,V);
eRef_g = real(D(i1,i1));

[V,D] = eig(NaD2_cpl.ops.H_NaD2_hyperfine);
[~,i2] = evec_ind({'f_NaD2','m_f_NaD2'},ref_e,NaD2_cpl,V);
eRef_e = real(D(i2,i2));

refDelta = eRef_e - eRef_g;

%% Find energy for a particular high field transition 
trs_g = [2,2]; %[F,mF] in coupled basis
trs_e = [1/2,3/2]; % [mJ,mI] in uncoupled basis

[V,D] = eig(NaGnd_cpl.ops.H_Na_hyperfine + NaGnd_cpl.ops.H0_Na_zeeman.*B);
[~,i1] = evec_ind({'f_Na','m_f_Na'},trs_g,NaGnd_cpl,V);
eTrs_g = real(D(i1,i1));

[VD2,DD2] = eig(NaD2.ops.H_NaD2_hyperfine + NaD2.ops.H0_NaD2_zeeman.*B);
[~,i2] = evec_ind({'m_j_NaD2','m_i_NaD2'},trs_e,NaD2,V);
allEs = diag(DD2);
eTrs_e = real(DD2(i2,i2));

transE = eTrs_e - eTrs_g;
deltaE = (transE - refDelta)/c.h/1e9;
disp([num2str(trs_g(1)),',',num2str(trs_g(2)),' -> ',num2str(trs_e(1)),',',num2str(trs_e(2)),' transition energy at B = ',num2str(866),'G is ',num2str(deltaE),' GHz above zero field 44->55'])

%% Plot transition strengths and energies for a particular starting state
trs_g = [2,2]; %[F,mF] in coupled basis
pol = [-1,0,1]; %Polarization, where 1 is sigma+, 0 is pi and -1 is sigma-

[V,D] = eig(NaGnd_cpl.ops.H_Na_hyperfine + NaGnd_cpl.ops.H0_Na_zeeman.*B);
[~,i1] = evec_ind({'f_Na','m_f_Na'},trs_g,NaGnd_cpl,V);
eTrs_g = real(D(i1,i1));
nonZeroInd = find(abs(V(:,i1))>0);
all_F = NaGnd_cpl.qnums.f_Na(nonZeroInd);
all_mF = NaGnd_cpl.qnums.m_f_Na(nonZeroInd);
Fpopns = abs(V(nonZeroInd,i1)).^2;

%Find energies relative to zero field transition for all eigenstates
transE = allEs - eTrs_g;
deltaE = (transE - refDelta)/c.h/1e9;
leg = {};

%Find all accessible excited states for each polarizatioin, plot energy vs
%transition strength
for j = 1:length(pol)
    allRelStr = zeros(size(NaD2.qnums,1),1);
    for i = 1:length(Fpopns)
        allRelStr = allRelStr + Fpopns(i)*relStr(all_F(i),all_mF(i),NaD2.qnums.j_NaD2(:),NaD2.qnums.m_j_NaD2(:),NaD2.qnums.i_NaD2(:),NaD2.qnums.m_i_NaD2(:),pol(j));
    end
    allEvalRelStr = abs(VD2).^2*allRelStr;
    figure(3)
    plt = scatter(deltaE,allEvalRelStr);
    set(plt,'markerFaceColor',plt.MarkerEdgeColor)
    set(plt,'markerFaceAlpha',0.5)
    xlabel('Energy [GHz] relative to 11->22 @ zero field')
    ylabel('Transition strength [arb]')
    hold on
    leg{end+1} = ['q = ',num2str(pol(j))];
end
legend(leg)
title(['F=',num2str(trs_g(1)),' m_F=',num2str(trs_g(2))])
hold off

%% Plot branching ratios for a particular upper state
trs_e = [1/2,3/2]; %[mJ,mI] in uncoupled basis
pol = [-1,0,1]; %Polarization, where 1 is sigma+, 0 is pi and -1 is sigma-

%Lists of mJ and mI with their admixtures in the eigenstate which is
%dominated by the mJ, mI given in trs_e
[~,i2] = evec_ind({'m_j_NaD2','m_i_NaD2'},trs_e,NaD2,VD2);
allEs = diag(DD2);
eTrs_e = real(DD2(i2,i2));
nonZeroInd = find(abs(VD2(:,i2))>0);
all_mJ = NaD2.qnums.m_j_NaD2(nonZeroInd);
all_mI = NaD2.qnums.m_i_NaD2(nonZeroInd);
uncpopns = abs(VD2(nonZeroInd,i2)).^2;

%Find transition strength of transitions to each ground eigenstate, avgd over
%polarizations
allEvalRelStr = zeros(size(NaGnd_cpl.qnums,1),1);
for j = 1:length(pol)
    allRelStr = zeros(size(NaGnd_cpl.qnums,1),1);
    for i = 1:length(uncpopns)
        allRelStr = allRelStr + uncpopns(i)*relStr(NaGnd_cpl.qnums.f_Na(:),NaGnd_cpl.qnums.m_f_Na(:),3/2,all_mJ(i),3/2,all_mI(i),pol(j));
    end
    allEvalRelStr = allEvalRelStr + 1/3*(abs(V).^2*allRelStr);
end

%Sort ground eigenstates back into order of qnums table, generate table of
%coupled ground eigenstates and their dominant quantum numbers
[~,V_ind] = max(abs(V(:,:,end)).^2,[],1);
allEvalRelStr = allEvalRelStr(V_ind);
[s,ind] = sort(allEvalRelStr,'descend');
brchTbl = NaGnd_cpl.qnums(ind(s>1e-6),3:4);
ratio = s(s>1e-6)/sum(s(s>1e-6));
brchTbl = addvars(brchTbl,ratio,'after','m_f_Na');
disp([num2str(trs_e(1)),',',num2str(trs_e(2)),' decay channels:'])
disp(brchTbl)

function res = relStr(F,mF,J,mJ,I,mI,q)
    J0 = 1/2;
    mJ0 = mF - mI;
    res = (2*F+1).*(w3j(J0,mJ0,I,mI,F,-mF).*w3j(J0,mJ0,1,q,J,-mJ)).^2;
end