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
qnums_ignore = {'S','Sigma','Lambda'};
fdata.qnums(:,qnums_ignore) = [];
cdata.qnums(:,qnums_ignore) = [];
Xdata.qnums(:,qnums_ignore) = [];

%% transition dipole moments
p = -1:1; % spherical index

% laser polarizations in cartesian coords
pol922 = [1 -1i 0]/sqrt(2); % sigma+
pol635 = [1 -1i 0]/sqrt(2); % sigma+

% laser spherical tensor operators
T922 = sphten(pol922);
T635 = sphten(pol635);

% transition dipole moments
TDM_c_f = operator_matrix(@transition_dipole_case_a,...
    {cdata.qnums,fdata.qnums},{'eta','J','Omega','m_J'},p);
TDM_X_c = operator_matrix(@transition_dipole_case_a,...
    {Xdata.qnums,cdata.qnums},{'eta','J','Omega','m_J'},p);

% off-diagonal Hamiltonian blocks
H_c_f = sum((-1).^reshape(p,1,1,3) .* TDM_c_f .* reshape(T922(p+2),1,1,[]),3);
H_X_c = sum((-1).^reshape(p,1,1,3) .* TDM_X_c .* reshape(T635(p+2),1,1,[]),3);

