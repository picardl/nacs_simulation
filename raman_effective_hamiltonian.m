clear;

c = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aUC';

power922 = 0.1; % W
power635 = 0.001; % W

waist = 13e-6; % m

% laser polarizations in cartesian coords
% sigma+
pol922 = [1 -1i 0]/sqrt(2);
pol635 = [1 -1i 0]/sqrt(2);

% sigma-
% pol922 = [1 1i 0]/sqrt(2);
% pol635 = [1 1i 0]/sqrt(2);

% % sigma+/-
% pol922 = [1 0 0];
% pol635 = [1 0 0];

% pi
% pol922 = [0 0 1];
% pol635 = [0 0 1];

%% laser stuff
% electric fields
E922 = sqrt(4*c.eta0*power922/(pi*waist^2)); % V/m
E635 = sqrt(4*c.eta0*power635/(pi*waist^2)); % V/m

% laser spherical tensor operators
T922 = sphten(pol922);
T635 = sphten(pol635);

%% load data
fdata = load(['feshbach_state_855G_' basis '.mat']);
cdata = load(['c3Sigma_state_855G_' basis '.mat']);
Xdata = load(['X1Sigma_state_855G_' basis '.mat']);

%% ignore quantum numbers that are not common to all 3 states
common_qnums = intersect(intersect(Xdata.qnums.Properties.VariableNames,...
    cdata.qnums.Properties.VariableNames),fdata.qnums.Properties.VariableNames);

fdata.qnums(:,~ismember(fdata.qnums.Properties.VariableNames,common_qnums)) = [];
cdata.qnums(:,~ismember(cdata.qnums.Properties.VariableNames,common_qnums)) = [];
Xdata.qnums(:,~ismember(Xdata.qnums.Properties.VariableNames,common_qnums)) = [];

% cut q nums we don't care about
qnums_ignore = {'Sigma','Lambda'};
fdata.qnums(:,qnums_ignore) = [];
cdata.qnums(:,qnums_ignore) = [];
Xdata.qnums(:,qnums_ignore) = [];


%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q = sphten([1 0 0]);

% transition dipole moments
switch basis
    case 'aFC'
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_aFC,...
            {cdata.qnums,fdata.qnums},{'eta','J','Omega','I','F','m_F'},p,T1q);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_aFC,...
            {Xdata.qnums,cdata.qnums},{'eta','J','Omega','I','F','m_F'},p,T1q);
    case 'aUC'
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
            {cdata.qnums,fdata.qnums},{'eta','J','Omega','m_J'},p,T1q);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_a,...
            {Xdata.qnums,cdata.qnums},{'eta','J','Omega','m_J'},p,T1q);
end

% off-diagonal Hamiltonian blocks
rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T922(p+2),1,1,3),3);
rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T635(p+2),1,1,3),3);

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

r = cdata.r;
psi_feshbach_interp = interp1(fdata.r,fdata.psi(:,:,1)',r,'spline')';
psi_X1Sigma_interp = interp1(Xdata.r,Xdata.psi_r,r,'spline');

a3S_c3S_elecTDM = (sqrt(1-c.c3Sigma.singlet_fraction)*(fdata.qnums.S==1) ...
    + sqrt(c.c3Sigma.singlet_fraction)*(fdata.qnums.S==0))...
    .*interp1(a3S_c3S_elecTDM_data(:,1)*c.abohr,a3S_c3S_elecTDM_data(:,2)*c.e*c.abohr,r);
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;

X1S_B1P_elecTDM = sqrt(c.c3Sigma.singlet_fraction).*...
    interp1(X1S_B1P_elecTDM_data(:,1)*c.abohr,X1S_B1P_elecTDM_data(:,2)*c.e*c.abohr,r);
X1S_B1P_elecTDM(isnan(X1S_B1P_elecTDM)) = 0;

vibronic_TDM_c_f = trapz(r,cdata.psi_r.*a3S_c3S_elecTDM.*psi_feshbach_interp,2);
vibronic_TDM_X_c = trapz(r,cdata.psi_r.*X1S_B1P_elecTDM.*psi_X1Sigma_interp,2);

H_TDM_c_f = rot_TDM_c_f.*vibronic_TDM_c_f' * E922;
H_TDM_X_c = rot_TDM_X_c.*vibronic_TDM_X_c' * E635;

%% let's plot some spectra
E_922 = 1e9*c.h*(325111 + linspace(5,18,1e3));

response_922 = sum(sum(abs(cdata.psi'*H_TDM_c_f).^2,2)./(cdata.E-fdata.E(1)-E_922 - 1i*c.c3Sigma.Gamma),1);

x = E_922*1e-9/c.h - 325111;
plot(x,imag(response_922),'linewidth',2)
xlim([min(x) max(x)])


cfnz = find(rot_TDM_c_f);
[r,c] = ind2sub(size(rot_TDM_c_f),cfnz);

switch basis
    case 'aUC'
        qnums_source = fdata.qnums(c,:);
        qnums_target = cdata.qnums(r,:);
        qnums_target.m_F = sum(qnums_target{:,{'m_i_Na','m_i_Cs','m_J'}},2);
        qnums_source
        qnums_target
    case 'aFC'
        qnums_source = fdata.qnums(c,:);
        qnums_target = cdata.qnums(r,:);
        qnums_source
        qnums_target
end

