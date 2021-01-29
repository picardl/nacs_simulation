clear;

const = constants();

% build an effective hamiltonian for the interaction between the feshbach
% molecule state, c3sigma, and x1sigma.

basis = 'aFC';

% power_up = 1; % W
% power_dn = 1; % W

% waist = 13e-6; % m

pol_up = 'sigp';
pol_dn = 'sigp';

%% laser stuff
% electric fields
Efield_up = 1; %sqrt(4*const.eta0*power_up/(pi*waist^2)); % V/m
Efield_dn = 1; %sqrt(4*const.eta0*power_dn/(pi*waist^2)); % V/m

% laser spherical tensor operators
switch pol_up
    case 'sigp'
        p_up = [1 -1i 0]/sqrt(2);
    case 'sigm'
        p_up = [1 1i 0]/sqrt(2);
    case 'pi'
        p_up = [0 0 1];
end
T_up = sphten(p_up);
switch pol_dn
    case 'sigp'
        p_dn = [1 -1i 0]/sqrt(2);
    case 'sigm'
        p_dn = [1 1i 0]/sqrt(2);
    case 'pi'
        p_dn = [0 0 1];
end
T_dn = sphten(p_dn);

%% load data
f = load(['data/feshbach_state_855G_' basis '.mat']);
c = load(['data/c3Sigma_state_855G_' basis '.mat']);
X = load(['data/X1Sigma_state_855G_' basis '.mat']);

%% ignore quantum numbers that are not common to all 3 states
common_qnums = intersect(intersect(X.qnums.Properties.VariableNames,...
    c.qnums.Properties.VariableNames),f.qnums.Properties.VariableNames);

f.qnums(:,~ismember(f.qnums.Properties.VariableNames,common_qnums)) = [];
c.qnums(:,~ismember(c.qnums.Properties.VariableNames,common_qnums)) = [];
X.qnums(:,~ismember(X.qnums.Properties.VariableNames,common_qnums)) = [];

% ignore some quantum numbers in rotational matrix element
qnums_ignore = {'S','Sigma','Lambda'};
f_ignore = ismember(f.qnums.Properties.VariableNames,qnums_ignore);
c_ignore = ismember(c.qnums.Properties.VariableNames,qnums_ignore);
X_ignore = ismember(X.qnums.Properties.VariableNames,qnums_ignore);

%% find energy reference point
switch basis
    case 'aFC'
        [~,ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],c,c.psi);
    case 'aUC'
        [~,ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],c,c.psi);
    case 'aIC'
        [~,ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
end
E0_up = mean(c.E(ind)-min(f.E));
E0_dn = mean(c.E(ind)-mean(X.E));

%% transition dipole moments, angular momentum part
p = -1:1; % spherical index

T1q_c_f = sphten([0 0 1]);
T1q_X_c = sphten([1 0 0]);

% transition dipole moments
switch basis
    case 'aFC'
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_aFC,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_c_f);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_aFC,...
            {c.qnums(:,~c_ignore),X.qnums(:,~X_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_X_c);
    case {'aUC','aIC'}
        rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
            {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_c_f);
        rot_TDM_X_c_components = operator_matrix(@transition_dipole_case_a,...
            {c.qnums(:,~c_ignore),X.qnums(:,~X_ignore)},{'eta','J','Omega','m_J'},p,T1q_X_c);
end

% off-diagonal Hamiltonian blocks
rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T_up(p+2),1,1,3),3);
rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T_dn(p+2),1,1,3),3);

%% transition dipole moments, vibronic part
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

a3S_c3S_elecTDM = (sqrt(1-const.c3Sigma.singlet_fraction)*(f.qnums.S==1)...
    .*interp1(a3S_c3S_elecTDM_data(:,1)*const.abohr,a3S_c3S_elecTDM_data(:,2)*const.e*const.abohr,c.r) ...
    + sqrt(const.c3Sigma.singlet_fraction)*(f.qnums.S==0)...
    .*interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,c.r));
a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;
rovib_TDM_c_f = rot_TDM_c_f.*permute(a3S_c3S_elecTDM,[3 1 2]);

X1S_c3S_elecTDM_interp = sqrt(const.c3Sigma.singlet_fraction)...
    *interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,X.r,'spline');

%% laser hamiltonian matrix elements
psi_c = c.psi.*permute(c.psi_r,[1 3 2]);
psi_feshbach_interp = permute(interp1(f.r,f.psi(:,:,1)',c.r,'spline')',[1 3 2]);
TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
TDM_c_f = trapz(c.r,TDM_R_c_f,3);
H_up = TDM_c_f*Efield_up;

psi_c_interp = interp1(c.r,c.psi_r,X.r,'spline');
vibronic_TDM_X_c = trapz(X.r,X.psi_r.*X1S_c3S_elecTDM_interp.*psi_c_interp);
H_dn = vibronic_TDM_X_c * Efield_dn * (c.psi'*rot_TDM_X_c*X.psi);

%% rabi frequencies per sqrt(mW)
% switch basis
%     case 'aFC'
%         [q_ind,v_ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],c,c.psi);
%     case 'aUC'
%         [q_ind,v_ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],c,c.psi);
%     case 'aIC'
%         [q_ind,v_ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
% end
% rabi_up_sqrt_mW = max(abs(H_up(v_ind)))/const.h*1e-6 *sqrt(1e-3/power_up);
% fprintf('up leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_up_sqrt_mW)
% 
% switch basis
%     case 'aFC'
%         [~,v_ind2] = evec_ind({'J','I','F','m_F'},[0,5,5,4],X,X.psi);
%     case 'aIC'
%         [~,v_ind2] = evec_ind({'J','m_J','I','m_I'},[0,0,5,4],X,X.psi);
%     case 'aUC'
%         [~,v_ind2] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[0,0,3/2,5/2],X,X.psi);
% end
% rabi_dn_sqrt_mW = max(abs(H_dn(v_ind,v_ind2)))/const.h*1e-6 *sqrt(1e-3/power_dn);
% fprintf('down leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_dn_sqrt_mW)

%% cut states that aren't coupled
% cut X states that aren't accessed via raman
X_cut = sqrt(abs(H_dn'*H_up))/const.h < 0.01*max(sqrt(abs(H_dn'*H_up))/const.h);
X.E(X_cut) = [];
X.psi(:,X_cut) = [];
H_dn(:,X_cut) = [];

% cut qnums in the ground state that don't have significant admixtures
X_qnums_cut = all(abs(X.psi).^2<1e-3,2);
X.qnums(X_qnums_cut,:) = [];
X.psi(X_qnums_cut,:) = [];
X.psi = X.psi./sqrt(sum(abs(X.psi).^2,1));

% cut c states that don't couple strongly to either ground state
c_cut = (all(abs(H_dn) < 0.1*max(abs(H_dn(:))),2)) & (abs(H_up) < 0.1*max(abs(H_up)));
c.psi(:,c_cut) = [];
c.E(c_cut) = [];
H_up(c_cut,:) = [];
H_dn(c_cut,:) = [];

% cut qnums in the excited state that don't have significant admixtures
c_qnums_cut = all(abs(c.psi).^2<1e-3,2);
c.qnums(c_qnums_cut,:) = [];
c.psi(c_qnums_cut,:) = [];
c.psi = c.psi./sqrt(sum(abs(c.psi).^2,1));

%% build effective hamiltonian
Nf = 1;
Nc = size(c.psi,2);
NX = size(X.psi,2);

Nstates = Nf + Nc + NX;

b1.qnums = build_basis({'eta','i'},{1,1},[0,0]);
b2.qnums = build_basis({'eta','i'},{2,1:size(c.psi,2)},[0,0]);
b3.qnums = build_basis({'eta','i'},{3,1:size(X.psi,2)},[0,0]);
b = join_basis(join_basis(b1,b2),b3);

b.ops.H0 = diag(cat(1,0,c.E-E0_up-f.E,X.E-(E0_up-E0_dn)-f.E));
b.ops.Hu = [zeros(Nf) H_up' zeros(Nf,NX); H_up zeros(Nc,Nc+NX); zeros(NX,Nstates)];
b.ops.Hd = [zeros(Nf,Nstates); zeros(Nc,Nf+Nc) H_dn; zeros(NX,Nf) H_dn' zeros(NX)];

b.ops.If = diag(cat(1,ones(Nf,1),zeros(Nc,1),zeros(NX,1)));
b.ops.Ic = diag(cat(1,zeros(Nf,1),ones(Nc,1),zeros(NX,1)));
b.ops.IX = diag(cat(1,zeros(Nf,1),zeros(Nc,1),ones(NX,1)));



%% save effective hamiltonian
fname = ['data/transfer_Heff_' basis '_' pol_up '_' pol_dn '.mat'];
save(fname,'b');
disp(['saved file ' fname])


