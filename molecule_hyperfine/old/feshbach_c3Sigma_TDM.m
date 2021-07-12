clear;

const = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aUC';

B = 855e-4;

%% load data
f = load(['../data/feshbach_state_' num2str(B*1e4) 'G_' basis '.mat']);
c = load('../data/c3Sigma_basis.mat');

%% ignore quantum numbers that are not common to all 3 states
common_qnums = intersect(c.qnums.Properties.VariableNames,f.qnums.Properties.VariableNames);

f.qnums(:,~ismember(f.qnums.Properties.VariableNames,common_qnums)) = [];
c.qnums(:,~ismember(c.qnums.Properties.VariableNames,common_qnums)) = [];

% ignore some quantum numbers in rotational matrix element
qnums_ignore = {'S','Sigma','Lambda'};
f_ignore = ismember(f.qnums.Properties.VariableNames,qnums_ignore);
c_ignore = ismember(c.qnums.Properties.VariableNames,qnums_ignore);

%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q_c_f = sphten([0 0 1]);

% transition dipole moments
switch basis
    case 'aFC'
        rot_TDM_c_f = operator_matrix(@transition_dipole_case_aFC,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_c_f);
    case {'aUC','aIC'}
        rot_TDM_c_f = operator_matrix(@transition_dipole_case_a,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_c_f);
end

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

a3S_c3S_elecTDM = (sqrt(1-const.c3Sigma.singlet_fraction)*(f.qnums.S==1)...
    .*interp1(a3S_c3S_elecTDM_data(:,1)*const.abohr,a3S_c3S_elecTDM_data(:,2)*const.e*const.abohr,c.r) ...
    + sqrt(const.c3Sigma.singlet_fraction)*(f.qnums.S==0)...
    .*interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,c.r));
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;
rovib_TDM_c_f = permute(rot_TDM_c_f,[1 2 4 3]).*permute(a3S_c3S_elecTDM,[3 1 2]);

%% laser hamiltonian matrix elements
psi_c = permute(c.psi_r,[1 3 2]);
psi_feshbach_interp = permute(interp1(f.r,f.psi(:,:,1)',c.r,'spline')',[1 3 2]);
TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
TDM = trapz(c.r,TDM_R_c_f,3);
TDM = permute(TDM,[1 4 2 3]);

%% save data
fname = ['../data/feshbach_c3Sigma_TDM_' num2str(B*1e4) 'G.mat'];
save(fname,'TDM');
