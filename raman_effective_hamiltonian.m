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

%% transition dipole moments

fdata;
