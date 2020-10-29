clear;

c = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aIC';

power_up = 50e-3; % W
power_dn = 80e-6; % W

waist = 13e-6; % m

% E0_up = 325110.965e9*c.h; 
% E0_dn = 472147*1e9*c.h; 

% laser polarizations in cartesian coords
% sigma+
pol_up = [1 -1i 0]/sqrt(2);
pol_dn = [1 -1i 0]/sqrt(2);

% sigma-
% pol_up = [1 1i 0]/sqrt(2);
% pol_dn = [1 1i 0]/sqrt(2);

% % sigma+/-
% pol_up = [1 0 0];
% pol_dn = [0 1 0];

% pi
% pol_up = [0 0 1];
% pol_dn = [0 0 1];

%% laser stuff
% electric fields
Efield_up = sqrt(4*c.eta0*power_up/(pi*waist^2)); % V/m
Efield_dn = sqrt(4*c.eta0*power_dn/(pi*waist^2)); % V/m

% laser spherical tensor operators
T_up = sphten(pol_up);
T_dn = sphten(pol_dn);

%% load data
fdata = load(['data/feshbach_state_855G_' basis '.mat']);
cdata = load(['data/c3Sigma_state_855G_' basis '.mat']);
Xdata = load(['data/X1Sigma_state_855G_' basis '.mat']);

%% ignore quantum numbers that are not common to all 3 states
common_qnums = intersect(intersect(Xdata.qnums.Properties.VariableNames,...
    cdata.qnums.Properties.VariableNames),fdata.qnums.Properties.VariableNames);

fdata.qnums(:,~ismember(fdata.qnums.Properties.VariableNames,common_qnums)) = [];
cdata.qnums(:,~ismember(cdata.qnums.Properties.VariableNames,common_qnums)) = [];
Xdata.qnums(:,~ismember(Xdata.qnums.Properties.VariableNames,common_qnums)) = [];

% ignore some quantum numbers in rotational matrix element
qnums_ignore = {'S','Sigma','Lambda'};
f_ignore = ismember(fdata.qnums.Properties.VariableNames,qnums_ignore);
c_ignore = ismember(cdata.qnums.Properties.VariableNames,qnums_ignore);
X_ignore = ismember(Xdata.qnums.Properties.VariableNames,qnums_ignore);

%% find energy reference point
switch basis
    case 'aFC'
        [~,ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],cdata,cdata.psi);
    case 'aUC'
        [~,ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],cdata,cdata.psi);
    case 'aIC'
        [~,ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],cdata,cdata.psi);
end
E0_up = mean(cdata.E(ind)-min(fdata.E));
E0_dn = mean(cdata.E(ind)-mean(Xdata.E));

%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q_c_f = sphten([0 0 1]);
T1q_X_c = sphten([1 0 0]);

% transition dipole moments
switch basis
    case 'aFC'
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_aFC,...
            {cdata.qnums(:,~c_ignore),fdata.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_c_f);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_aFC,...
            {cdata.qnums(:,~c_ignore),Xdata.qnums(:,~X_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_X_c);
    case {'aUC','aIC'}
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
            {cdata.qnums(:,~c_ignore),fdata.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_c_f);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_a,...
            {cdata.qnums(:,~c_ignore),Xdata.qnums(:,~X_ignore)},{'eta','J','Omega','m_J'},p,T1q_X_c);
end

% off-diagonal Hamiltonian blocks
rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T_up(p+2),1,1,3),3);
rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T_dn(p+2),1,1,3),3);

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

a3S_c3S_elecTDM = (sqrt(1-c.c3Sigma.singlet_fraction)*(fdata.qnums.S==1) ...
    + sqrt(c.c3Sigma.singlet_fraction)*(fdata.qnums.S==0))...
    .*interp1(a3S_c3S_elecTDM_data(:,1)*c.abohr,a3S_c3S_elecTDM_data(:,2)*c.e*c.abohr,cdata.r);
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;
rovib_TDM_c_f = rot_TDM_c_f.*permute(a3S_c3S_elecTDM,[3 1 2]);

X1S_c3S_elecTDM_interp = interp1(X1S_B1P_elecTDM_data(:,1)*c.abohr,X1S_B1P_elecTDM_data(:,2)*c.e*c.abohr,Xdata.r,'spline');

%% laser hamiltonian matrix elements
psi_c = cdata.psi.*permute(cdata.psi_r,[1 3 2]);
psi_feshbach_interp = permute(interp1(fdata.r,fdata.psi(:,:,1)',cdata.r,'spline')',[1 3 2]);
TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
TDM_c_f = trapz(cdata.r,TDM_R_c_f,3);
H_up = TDM_c_f*Efield_up;

psi_c_interp = interp1(cdata.r,cdata.psi_r,Xdata.r,'spline');
vibronic_TDM_X_c = sqrt(c.c3Sigma.singlet_fraction)*trapz(Xdata.r,Xdata.psi_r.*X1S_c3S_elecTDM_interp.*psi_c_interp);
H_dn = vibronic_TDM_X_c * Efield_dn * (cdata.psi'*rot_TDM_X_c*Xdata.psi);

%% plot scattering rates
x = linspace(-10,10,1e3);
E_up = E0_up + 1e9*c.h*x;
E_dn = E0_dn + 1e9*c.h*x;

decay_up = (c.c3Sigma.Gamma/(2*c.h)) * sum(2*abs(H_up).^2./(2*abs(H_up).^2 + 4*(cdata.E - fdata.E' - E_up).^2 + c.c3Sigma.Gamma^2),1);
decay_dn = (c.c3Sigma.Gamma/(2*c.h)) * sum(2*sum(abs(H_dn).^2,2)./(2*sum(abs(H_dn).^2,2) + 4*(cdata.E - mean(Xdata.E) - E_dn).^2 + c.c3Sigma.Gamma^2),1);

figure(1);
clf;
hold on;
box on;
plot(x,decay_up*1e-6,'linewidth',2);
plot(x,decay_dn*1e-6,'linewidth',2);
hold off;
set(gca,'fontsize',14)
xlabel('detuning (GHz)')
ylabel('scattering rate (MHz)');
xlim([min(x) max(x)])
% ylim([0 0.01])

%% rabi frequencies per sqrt(mW)
switch basis
    case 'aFC'
        [q_ind,v_ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],cdata,cdata.psi);
    case 'aUC'
        [q_ind,v_ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],cdata,cdata.psi);
    case 'aIC'
        [q_ind,v_ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],cdata,cdata.psi);
end
rabi_up_sqrt_mW = max(abs(H_up(v_ind)))/c.h*1e-6 *sqrt(1e-3/power_up);
fprintf('up leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_up_sqrt_mW)

switch basis
    case 'aFC'
        [~,v_ind2] = evec_ind({'J','I','F','m_F'},[0,5,5,4],Xdata,Xdata.psi);
    case 'aIC'
        [~,v_ind2] = evec_ind({'J','m_J','I','m_I'},[0,0,5,4],Xdata,Xdata.psi);
    case 'aUC'
        [~,v_ind2] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[0,0,3/2,5/2],Xdata,Xdata.psi);
end
rabi_dn_sqrt_mW = max(abs(H_dn(v_ind,v_ind2)))/c.h*1e-6 *sqrt(1e-3/power_dn);
fprintf('down leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_dn_sqrt_mW)

%% cut states
% cut X states that aren't accessed via raman
X_cut = sqrt(abs(H_dn'*H_up))/c.h*1e-6 < 0.1;
Xdata.E(X_cut) = [];
Xdata.psi(:,X_cut) = [];
H_dn(:,X_cut) = [];

% cut c states that don't couple strongly to either ground state
c_cut = (all(abs(H_dn) < 0.1*max(abs(H_dn(:))),2)) & (abs(H_up) < 0.1*max(abs(H_up)));
cdata.psi(:,c_cut) = [];
cdata.E(c_cut) = [];
H_up(c_cut,:) = [];
H_dn(c_cut,:) = [];

%% build effective hamiltonian
Nf = 1;
Nc = size(cdata.psi,2);
NX = size(Xdata.psi,2);

Nstates = Nf + Nc + NX;

H0 = diag(cat(1,0,cdata.E-E0_up-fdata.E-1i*c.c3Sigma.Gamma,Xdata.E-(E0_up-E0_dn)-fdata.E));
H_up = [zeros(Nf) H_up' zeros(Nf,NX); H_up zeros(Nc,Nc+NX); zeros(NX,Nstates)];
H_dn = [zeros(Nf,Nstates); zeros(Nc,Nf+Nc) H_dn; zeros(NX,Nf) H_dn' zeros(NX)];

(H0+H_up+H_dn)/c.h * 1e-9

sqrt(abs(H_dn'*H_up))/c.h*1e-6
