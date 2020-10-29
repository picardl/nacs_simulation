clear;

c = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aUC';

power_up = 100e-3; % W
power_dn = 60e-6; % W

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
% pol_dn = [1 0 0];

% pi
% pol_up = [0 0 1];
% pol_dn = [0 0 1];

%% laser stuff
% electric fields
E_up = sqrt(4*c.eta0*power_up/(pi*waist^2)); % V/m
E_dn = sqrt(4*c.eta0*power_dn/(pi*waist^2)); % V/m

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

% cut q nums we don't care about
qnums_ignore = {'S','Sigma','Lambda'};
f_ignore = ismember(fdata.qnums.Properties.VariableNames,qnums_ignore);
c_ignore = ismember(cdata.qnums.Properties.VariableNames,qnums_ignore);
X_ignore = ismember(Xdata.qnums.Properties.VariableNames,qnums_ignore);

%% find energy of nuclear spin stretched m_J=0 state
[~,ind] = max(evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,0,3/2,7/2],cdata,cdata.psi)==evec_leading_percentages(cdata.psi,1),[],2);
E0_up = mean(cdata.E(ind)-min(fdata.E));
E0_dn = mean(cdata.E(ind)-mean(Xdata.E));

%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q_c_f = sphten([1 0 0]);
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

X1S_c3S_elecTDM = sqrt(c.c3Sigma.singlet_fraction)*mean(interp1(X1S_B1P_elecTDM_data(:,1)*c.abohr,X1S_B1P_elecTDM_data(:,2)*c.e*c.abohr,Xdata.r,'spline'));

%% laser effective hamiltonian matrix elements
psi_c = cdata.psi.*permute(cdata.psi_r,[1 3 2]);

psi_feshbach_interp = permute(interp1(fdata.r,fdata.psi(:,:,1)',cdata.r,'spline')',[1 3 2]);
TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
TDM_c_f = trapz(cdata.r,TDM_R_c_f,3);
H_up = TDM_c_f*E_up;

psi_c_interp = interp1(cdata.r,cdata.psi_r,Xdata.r,'spline');
vib_TDM_X_c = trapz(Xdata.r,Xdata.psi_r.*psi_c_interp);
H_dn = vib_TDM_X_c * X1S_c3S_elecTDM * E_dn * (cdata.psi'*rot_TDM_X_c*Xdata.psi);

%% let's plot some spectra
% E_up = E0_up + 1e9*c.h*linspace(-5,15,1e3);
% E_dn = E0_dn + 1e9*c.h*linspace(-5,15,1e3);

E_up = E0_up + 1e9*c.h*linspace(-25,-2,1e3);
E_dn = E0_dn + 1e9*c.h*linspace(-25,-2,1e3);

decay_up = (c.c3Sigma.Gamma/(2*c.h)) * sum(2*abs(H_up).^2./(2*abs(H_up).^2 + 4*(cdata.E - fdata.E' - E_up).^2 + c.c3Sigma.Gamma^2),1);
decay_dn = (c.c3Sigma.Gamma/(2*c.h)) * sum(2*sum(abs(H_dn).^2,2)./(2*sum(abs(H_dn).^2,2) + 4*(cdata.E - mean(Xdata.E) - E_dn).^2 + c.c3Sigma.Gamma^2),1);

figure(1);
clf;
hold on;
box on;
x1 = (E_up-E0_up)/c.h*1e-9;
plot(x1,decay_up*1e-6,'linewidth',2);
x2 = (E_dn-E0_dn)/c.h*1e-9;
plot(x2,decay_dn*1e-6,'linewidth',2);
hold off;
set(gca,'fontsize',14)
xlabel('detuning (GHz)')
ylabel('scattering rate (MHz)');
xlim([min(x1) max(x1)])


sqrt(abs(H_dn'*H_up))/c.h * 1e-6

%% build effective hamiltonian

% Nexcited_keep = 1;
% 
% [val_up,order_up] = sort(abs(H_up),'descend');
% 
% % Xraman = find(H_dn*H_up);
% % 
% % nz_dn = find(abs(H_dn(:)*1e-6/c.h));
% % nz_up = find(abs(H_up(:)*1e-6/c.h));
% % 
% % [val,order] = sort(abs(H_up(:)*1e-6/c.h),'descend');
% 
% h_dn_nz = find(abs(H_dn));
% [r_dn_nz,c_dn_nz] = ind2sub(size(H_dn),h_dn_nz);
% H_dn_nz = diag(H_dn(r_dn_nz,c_dn_nz));
% 
% [val_dn,order_dn] = sort(abs(H_dn_nz),'descend');
