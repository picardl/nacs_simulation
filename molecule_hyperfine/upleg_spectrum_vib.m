%Generate up-leg spectrum of a3Sigma - c3Sigma for a particular v'

const = constants;

B = (860)*1e-4; %B Field in Gauss. updated FR in model is at 864.1
initState = 1; %1 for FB mol, 2 for atoms
basis = 'aFC'; %Use fully coupled basis (F,mF,J,mJ)
recompute = 2;
Jmax = 2;%Max J to inlcude in rotational basis
mtot = [3:5]; %Range of mF to include in hyperfine basis
f_vib = 320000e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
pol_up = 'sigp';
pol_dn = 'sigp';%Polarization of up-leg, can be sigp, sigm, pi, or sigpm
% 
% %Compute effective Hamiltonian for up-leg (and down-leg) including
% %rovibrational contribution to TDM
raman_data = raman_effective_hamiltonian(B,basis,recompute,Jmax,mtot,f_vib);

% initial state
psi_init = zeros(size(raman_data.f.E));
%Select lowest energy state fromfeshbach hamiltonian, corresponding to the
%FB molecule
if initState == 1
    indState = find(raman_data.f.E < raman_data.f.E_lowest_chan_threshold);
    psi_init(indState(end)) = 1; %Find highest energy state below threshold
elseif initState == 2
    indState = find(raman_data.f.E > raman_data.f.E_lowest_chan_threshold);
    psi_init(indState(1)) = 1; %Find lowest energy state above threshold
end

% laser stuff for figuring out realistic scattering rates
power = 100e-3; %W
waist = 50e-6; %m
p_up = sphten(polSwitch(pol_up));
p_down = sphten(polSwitch(pol_dn));
t = 10e-3;

Efield = sqrt(4*const.eta0*power/(pi*waist^2));
excl = raman_data.c.qnums.Omega == 0;
Ndet = 1e4;
freq_offs = min(min(raman_data.c.E - raman_data.f.E))/const.h;
detuning = linspace(-1e9,13e9,Ndet);

%% 
%Calculate up-leg rabi rates to each state and corresponding total
%scattering rates
Gamma = const.h*15e6;
rabi = Efield*sum(raman_data.H_up.*psi_init.*permute(p_up,[1 3 4 2]),4);
delta = raman_data.c.E - raman_data.f.E - const.h*reshape(freq_offs + detuning,1,1,[]);
Rsc = permute(sum(sum((Gamma/(2*const.h)) * 2*abs(rabi).^2./(2*abs(rabi).^2 + 4*delta.^2 + Gamma.^2),2),1),[2 3 1]);
t = 3e-6;
popn = exp(-Rsc*t);

figure(1); clf;
subplot(3,1,1)
x = (raman_data.c.E - raman_data.f.E(initState))*1e-9/const.h - freq_offs*1e-9;
% x = x - min(min(raman_data.f.E/const.h/1e9));
y = abs(rabi(:,initState)).^2./max(abs(rabi(:,initState)).^2);
scatter(x,y)
xlabel('Up-leg transition freq [GHz]')
ylabel('Normalized transition strength')
xlim([-1,13])
hold on
line([x x]', [y zeros(size(y))]','Color',[0 0.4470 0.7410],'LineWidth',2)
title('v = 22, \sigma_{+-} polarization, B = FR - 4G')

subplot(3,1,2)
plot(detuning*1e-9,Rsc*1e-6)
xlabel('Up-leg transition freq [GHz]')
xlim([-1,13])
ylabel('R_{sc} (10^6 s^-1)')

subplot(3,1,3)
plot(detuning*1e-9,popn)
xlabel('Up-leg transition freq [GHz]')
xlim([-1,13])
ylabel('Atom popn')


%% Extract states
dispTbl = table();
dispTbl.J = round(0.5*(sqrt(4*diag(raman_data.c.psi'*raman_data.c.ops.J_sq*raman_data.c.psi)+1)-1));
dispTbl.mJ = round(diag(raman_data.c.psi'*raman_data.c.ops.J_z*raman_data.c.psi));
dispTbl.F = round(sqrt(diag(raman_data.c.psi'*raman_data.c.ops.F_sq*raman_data.c.psi)));
dispTbl.mF = round(diag(raman_data.c.psi'*raman_data.c.ops.F_z*raman_data.c.psi));
dispTbl.delta = x;
dispTbl.normRabi = y;
mS = diag(raman_data.c.psi'*raman_data.c.ops.S_z*raman_data.c.psi);
[srtRabi,ind] = sortrows(dispTbl,'normRabi','descend');
[~,maxBase] = max(raman_data.c.psi.^2,[],1);
srtc = raman_data.c.qnums(maxBase(ind),:);
% srtRabi(srtRabi.J == 3,:) = [];
disp(['States sorted by upleg transition strength for pol_up = ',pol_up])
disp(srtRabi(1:50,:))

if 0
    %%
    disp('The eigenstate with (rounded) quantum numbers')
    sortedRow = 6;
    disp(srtRabi(sortedRow,:))
    disp('Is composed of the following states in c3Sigma basis:')
    initRow = ind(sortedRow);
    psiRow = raman_data.c.psi(:,initRow);
    [sortPsi,indSP] = sort(psiRow.*conj(psiRow),'descend');
    componentQnums = raman_data.c.qnums(indSP,:);
    componentQnums.amp = sortPsi;
    componentQnums.psi = psiRow(indSP);
    disp(componentQnums(1:10,:))
end

nonRTable = table();
nonRTable.J = (0.5*(sqrt(4*diag(raman_data.c.psi'*raman_data.c.ops.J_sq*raman_data.c.psi)+1)-1));
nonRTable.mJ = (diag(raman_data.c.psi'*raman_data.c.ops.J_z*raman_data.c.psi));
nonRTable.F = (sqrt(diag(raman_data.c.psi'*raman_data.c.ops.F_sq*raman_data.c.psi)));
nonRTable.mF = (diag(raman_data.c.psi'*raman_data.c.ops.F_z*raman_data.c.psi));
nonRTable.I = (sqrt(diag(raman_data.c.psi'*raman_data.c.ops.I_sq*raman_data.c.psi)));
nonRTable.mI = (diag(raman_data.c.psi'*raman_data.c.ops.I_z*raman_data.c.psi));
nonRTable.mI_Na = (diag(raman_data.c.psi'*raman_data.c.ops.i_Na_z*raman_data.c.psi));
nonRTable.mI_Cs = (diag(raman_data.c.psi'*raman_data.c.ops.i_Cs_z*raman_data.c.psi));
nonRTable.delta = x;
nonRTable = nonRTable(ind,:);

%Down leg
rabiDn = Efield*sum(raman_data.H_dn.*reshape(p_down,1,1,1,1,3),5)/const.h;
% rabiDn = abs(rabiDn)/max(max(abs(rabiDn)));
rabiUp = rabi/const.h;%abs(rabi)/max(max(abs(rabi)));
ramanAmpProd = conj(rabiUp(:,initState)).*rabiDn;
ramanAmpAbs = abs(ramanAmpProd)/1e6/1e9;
upLegAmp = sum(ramanAmpAbs,2);
dnLegAmp = sum(ramanAmpAbs,1);
dispTbl.OmegaRamanMHzperGHzDet = upLegAmp;
dispTbl.rabiUp = abs(rabiUp(:,initState))/1e6;
dispTbl.crabiUp = rabiUp(:,initState)/1e6;
srtRamanUp = sortrows(dispTbl,'OmegaRamanMHzperGHzDet','descend');
srtRamanUp(srtRamanUp.J == 3,:) = [];
srtRamanUp.normRabiUp = srtRamanUp.rabiUp/(srtRamanUp.rabiUp(3));
% srtRamanUp.delta = srtRamanUp.delta - srtRamanUp.delta(3);
disp(['States sorted by raman transition strength summed over down leg for pol_up = ',pol_up,' pol_down = ',pol_dn])
disp(srtRamanUp(1:12,:))

if 0
    mulitCpl = srtRamanUp(:,[1:5,8:9]);
    NRows = size(mulitCpl,1);
    mulitCpl.rabiDn1272 = zeros(NRows,1);
    mulitCpl.rabiDn3252 = zeros(NRows,1);
    mulitCpl.absrabiDn1272 = zeros(NRows,1);
    mulitCpl.absrabiDn3252 = zeros(NRows,1);
    for i = 1:size(mulitCpl,1)
        indUp = find(sum(table2array(dispTbl(:,1:5))==table2array(srtRamanUp(i,1:5)),2)==5);
        dispDn = table();
        EGnd = diag(raman_data.X.psi'*raman_data.X.ops.H*raman_data.X.psi);
        dispDn.J = round(0.5*(sqrt(4*diag(raman_data.X.psi'*raman_data.X.ops.J_sq*raman_data.X.psi)+1)-1));
        dispDn.mJ = round(diag(raman_data.X.psi'*raman_data.X.ops.J_z*raman_data.X.psi));
        dispDn.F = round(sqrt(diag(raman_data.X.psi'*raman_data.X.ops.F_sq*raman_data.X.psi)));
        dispDn.mF = round(diag(raman_data.X.psi'*raman_data.X.ops.F_z*raman_data.X.psi));
        dispDn.mNa = round(2*diag(raman_data.X.psi'*raman_data.X.ops.F1_z*raman_data.X.psi))/2;
        dispDn.mCs = round(2*diag(raman_data.X.psi'*raman_data.X.ops.F2_z*raman_data.X.psi))/2;
        dispDn.EMHz = (EGnd - min(EGnd))/const.h/1e6;
        dispDn.OmegaRamanMHz = ramanAmpAbs(indUp,:)';
        dispDn.crabiDnMHz = rabiDn(indUp,:)'/1e6;
        dispDn.rabiDnMHz = abs(rabiDn(indUp,:)')/1e6;
        dispDn.ramanAmpProd = ramanAmpProd(indUp,:)';
        ind1272 = find(sum(table2array(dispDn(:,[1,5,6]))== [0,0.5,3.5],2) == 3);
        ind3253 = find(sum(table2array(dispDn(:,[1,5,6]))== [0,1.5,2.5],2) == 3);
        mulitCpl.rabiDn1272(i) = dispDn.crabiDnMHz(ind1272);
        mulitCpl.rabiDn3252(i) = dispDn.crabiDnMHz(ind3253);
        mulitCpl.absrabiDn1272(i) = abs(dispDn.crabiDnMHz(ind1272));
        mulitCpl.absrabiDn3252(i) = abs(dispDn.crabiDnMHz(ind3253));
    end

    % disp(mulitCpl(1:10,:))
    % writetable(mulitCpl,['../data/',datestr(now,'yyyymmdd'),'ramanCouplings_hf_up-',pol_up,'_dn-',pol_dn,'.csv']);
end
if 0
    %%
    indSortUp = 3; %Select which row from sorted raman strength to use as intermediate
    indUp = find(sum(table2array(dispTbl(:,1:5))==table2array(srtRamanUp(indSortUp,1:5)),2)==5);
    dispDn = table();
    EGnd = diag(raman_data.X.psi'*raman_data.X.ops.H*raman_data.X.psi);
    dispDn.J = round(0.5*(sqrt(4*diag(raman_data.X.psi'*raman_data.X.ops.J_sq*raman_data.X.psi)+1)-1));
    dispDn.mJ = round(diag(raman_data.X.psi'*raman_data.X.ops.J_z*raman_data.X.psi));
    dispDn.F = round(sqrt(diag(raman_data.X.psi'*raman_data.X.ops.F_sq*raman_data.X.psi)));
    dispDn.mF = round(diag(raman_data.X.psi'*raman_data.X.ops.F_z*raman_data.X.psi));
    dispDn.mNa = round(2*diag(raman_data.X.psi'*raman_data.X.ops.F1_z*raman_data.X.psi))/2;
    dispDn.mCs = round(2*diag(raman_data.X.psi'*raman_data.X.ops.F2_z*raman_data.X.psi))/2;
    dispDn.EMHz = (EGnd - min(EGnd))/const.h/1e6;
    dispDn.OmegaRamanMHz = ramanAmpAbs(indUp,:)';
    dispDn.rabiDnMHz = abs(rabiDn(indUp,:)')/1e6;
    dispDn.ramanAmpProd = ramanAmpProd(indUp,:)';
    dispDnSrt = sortrows(dispDn,'OmegaRamanMHz','descend');
    dispDnSrt(dispDnSrt.J ~= 0,:) = [];
    
    disp('With intermediate state:')
    disp(srtRamanUp(indSortUp,:))
    disp(['Rovib ground states sorted by raman transition strength for pol_up = ',pol_up,' pol_down = ',pol_dn])
    disp(dispDnSrt(1:6,:))
end
function pOut = polSwitch(pol)
    switch pol
        case 'sigp'
            pOut = [1 -1i 0]/sqrt(2);
        case 'sigm'
            pOut = [1 1i 0]/sqrt(2);
        case 'pi'
            pOut = [0 0 1];
        case 'sigpm'
            pOut = [1 0 0];
    end
end