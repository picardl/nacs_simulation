clear;

c = constants();

%% load data
fdata = load('feshbach_state_855G.mat');
cdata = load('c3Sigma_state_855G.mat');
Xdata = load('X1Sigma_state_855G.mat');

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

% laser polarizations in cartesian coords
pol922 = [1 -1i 0]/sqrt(2); % sigma+
pol635 = [1 -1i 0]/sqrt(2); % sigma+

% laser spherical tensor operators
T922 = sphten(pol922);
T635 = sphten(pol635);

% transition dipole moments
rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
    {cdata.qnums,fdata.qnums},{'eta','J','Omega','m_J'},p);
rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_a,...
    {Xdata.qnums,cdata.qnums},{'eta','J','Omega','m_J'},p);

% off-diagonal Hamiltonian blocks
rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T922(p+2),1,1,3),3);
rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T635(p+2),1,1,3),3);

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

r = cdata.r;
for i = 1:size(fdata.psi,3)
    psi_lower_feshbach(:,:,i) = interp1(fdata.r,fdata.psi(:,:,i)',r,'spline')';
end
psi_lower_X1Sigma = interp1(Xdata.r,Xdata.psi_r,r,'spline');

a3S_c3S_elecTDM = (sqrt(1-c.c3Sigma.singlet_fraction)*(fdata.qnums.S==1) + sqrt(c.c3Sigma.singlet_fraction)*(fdata.qnums.S==0)).*interp1(a3S_c3S_elecTDM_data(:,1)*c.abohr,a3S_c3S_elecTDM_data(:,2)*c.e*c.abohr,r);
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;

X1S_B1P_elecTDM = sqrt(c.c3Sigma.singlet_fraction).*interp1(X1S_B1P_elecTDM_data(:,1)*c.abohr,X1S_B1P_elecTDM_data(:,2)*c.e*c.abohr,r);
X1S_B1P_elecTDM(isnan(X1S_B1P_elecTDM)) = 0;

vibronic_TDM_c_f = trapz(r,cdata.psi_r.*a3S_c3S_elecTDM.*psi_lower_feshbach,2);
vibronic_TDM_X_c = trapz(r,cdata.psi_r.*X1S_B1P_elecTDM.*psi_lower_X1Sigma,2);

H_TDM_c_f = rot_TDM_c_f.*vibronic_TDM_c_f';
H_TDM_X_c = rot_TDM_X_c.*vibronic_TDM_X_c';

%% let's plot some spectra



