clear;

c = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aFC';

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
            {Xdata.qnums(:,~X_ignore),cdata.qnums(:,~c_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_X_c);
    case {'aUC','aIC'}
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
            {cdata.qnums(:,~c_ignore),fdata.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_c_f);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_a,...
            {Xdata.qnums(:,~X_ignore),cdata.qnums(:,~c_ignore)},{'eta','J','Omega','m_J'},p,T1q_X_c);
end

% off-diagonal Hamiltonian blocks
rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T922(p+2),1,1,3),3);
rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T635(p+2),1,1,3),3);

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

a3S_c3S_elecTDM = (sqrt(1-c.c3Sigma.singlet_fraction)*(fdata.qnums.S==1) ...
    + sqrt(c.c3Sigma.singlet_fraction)*(fdata.qnums.S==0))...
    .*interp1(a3S_c3S_elecTDM_data(:,1)*c.abohr,a3S_c3S_elecTDM_data(:,2)*c.e*c.abohr,cdata.r);
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;
rovib_TDM_c_f = rot_TDM_c_f.*permute(a3S_c3S_elecTDM,[3 1 2]);

X1S_c3S_elecTDM = mean(interp1(X1S_B1P_elecTDM_data(:,1)*c.abohr,X1S_B1P_elecTDM_data(:,2)*c.e*c.abohr,Xdata.r,'spline'));

%% laser effective hamiltonian matrix elements
psi_c = cdata.psi.*permute(cdata.psi_r,[1 3 2]);

psi_feshbach_interp = permute(interp1(fdata.r,fdata.psi(:,:,1)',cdata.r,'spline')',[1 3 2]);
TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
TDM_c_f = trapz(cdata.r,TDM_R_c_f,3);
H_922 = TDM_c_f*E922;

psi_c_interp = interp1(cdata.r,cdata.psi_r,Xdata.r,'spline');
vib_TDM_X_c = trapz(Xdata.r,Xdata.psi_r.*psi_c_interp);

H_635 = vib_TDM_X_c*X1S_c3S_elecTDM*E635 * (Xdata.psi'*rot_TDM_X_c*cdata.psi);

H_635'

%% let's plot some spectra
E_922 = 1e9*c.h*(325111 + linspace(5,18,1e3));

decay_rate = sum(2*abs(H_922).^2./(2*abs(H_922).^2 + 4*(cdata.E - fdata.E' - E_922).^2 + c.c3Sigma.Gamma^2),1);

plot(E_922/c.h*1e-9 - 325000,decay_rate);
set(gca,'fontsize',14)
xlabel('up leg detuning - 325000 (GHz)')
ylabel('scattering rate (MHz)');
title(['P_{922nm} = ' num2str(power922*1e3) ' mW'])

%%

