%Generate up-leg spectrum of a3Sigma - c3Sigma for a particular v'

const = constants;

B = (852.5 - 0)*1e-4; %B Field in Gauss. FR in model is at 852.5 G
initState = 2; %1 for FB mol, 2 for atoms
basis = 'aFC'; %Use fully coupled basis (F,mF,J,mJ)
recompute = 2;
Jmax = 3;%Max J to inlcude in rotational basis
mtot = [-6:6]; %Range of mF to include in hyperfine basis
f_vib = 320000e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
pol_up = 'sigpm';
pol_dn = 'sigpm';%Polarization of up-leg, can be sigp, sigm, pi, or sigpm
% 
% %Compute effective Hamiltonian for up-leg (and down-leg) including
% %rovibrational contribution to TDM
% raman_data = raman_effective_hamiltonian(B,basis,recompute,Jmax,mtot,f_vib);

if length(raman_data.f.E) > 2
    initState = initState + 1;
end

% initial state
%Select lowest energy state fromfeshbach hamiltonian, corresponding to the
%FB molecule
psi_init = zeros(size(raman_data.f.E));
psi_init(initState) = 1;

% laser stuff for figuring out realistic scattering rates
power = 100e-3; %W
waist = 50e-6; %m
p_up = sphten(polSwitch(pol_up));
p_down = sphten(polSwitch(pol_dn));
t = 10e-3;

Efield = sqrt(4*const.eta0*power/(pi*waist^2));

Ndet = 1e4;
freq_offs = (320064)*1e9;
detuning = linspace(0,12e9,Ndet);

%%
%Calculate up-leg rabi rates to each state and corresponding total
%scattering rates
Gamma = const.c3Sigma.Gamma/5;
rabi = Efield*sum(raman_data.H_up.*psi_init.*permute(p_up,[1 3 4 2]),4);
delta = raman_data.c.E - raman_data.f.E - const.h*reshape(freq_offs + detuning,1,1,[]);
Rsc = permute(sum(sum((Gamma/(2*const.h)) * 2*abs(rabi).^2./(2*abs(rabi).^2 + 4*delta.^2 + Gamma.^2),2),1),[2 3 1]);

figure(1); clf;
subplot(2,1,1)
x = raman_data.c.E*1e-9/const.h - freq_offs*1e-9;
x = x - min(min(raman_data.f.E/const.h/1e9));
y = abs(rabi(:,initState)).^2./max(abs(rabi(:,initState)).^2);
scatter(x,y)
xlabel('Up-leg transition freq [GHz]')
ylabel('Normalized transition strength')
xlim([0,12])
hold on
line([x x]', [y zeros(size(y))]','Color',[0 0.4470 0.7410],'LineWidth',2)
title('v = 22, \sigma_- polarization, B = FR - 6G')

subplot(2,1,2)
plot(detuning*1e-9,Rsc*1e-6)
xlabel('Up-leg transition freq [GHz]')
ylabel('R_{sc} (10^6 s^-1)')

%%Extract states
dispTbl = table();
dispTbl.J = round(0.5*(sqrt(4*diag(raman_data.c.psi'*raman_data.c.ops.J_sq*raman_data.c.psi)+1)-1));
dispTbl.mJ = round(diag(raman_data.c.psi'*raman_data.c.ops.J_z*raman_data.c.psi));
dispTbl.F = round(sqrt(diag(raman_data.c.psi'*raman_data.c.ops.F_sq*raman_data.c.psi)));
dispTbl.mF = round(diag(raman_data.c.psi'*raman_data.c.ops.F_z*raman_data.c.psi));
dispTbl.delta = x;
dispTbl.normRabi = y;
mS = diag(raman_data.c.psi'*raman_data.c.ops.S_z*raman_data.c.psi);
srtRabi = sortrows(dispTbl,'normRabi','descend');
srtRabi(srtRabi.J == 3,:) = [];
disp(['States sorted by upleg transition strength for pol_up = ',pol_up])
disp(srtRabi(1:8,:))

%Down leg
rabiDn = Efield*sum(raman_data.H_dn.*reshape(p_down,1,1,1,1,3),5)/const.h;
% rabiDn = abs(rabiDn)/max(max(abs(rabiDn)));
rabiUp = rabi/const.h;%abs(rabi)/max(max(abs(rabi)));
ramanAmpProd = conj(rabiUp(:,initState)).*rabiDn;
ramanAmpAbs = sqrt(abs(ramanAmpProd))/1e6;
upLegAmp = sum(ramanAmpAbs,2);
dnLegAmp = sum(ramanAmpAbs,1);
dispTbl.OmegaRamanSumUpMHz = upLegAmp;
srtRamanUp = sortrows(dispTbl,'OmegaRamanSumUpMHz','descend');
srtRamanUp(srtRamanUp.J == 3,:) = [];
disp(['States sorted by raman transition strength summed over down leg for pol_up = ',pol_up,' pol_down = ',pol_dn])
disp(srtRamanUp(1:8,:))

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