clear;

const = constants();

% Determine tdm dependence on B field close to FB resonance

basis = 'aFC';
%% load data
BFields = [852:0.5:858];
tdmUp = zeros(size(BFields));
tdmUpSum = zeros(size(BFields));
fAll = {};
cAll = {};
for i = 1:length(BFields)
    fAll{i} = load(['data/feshbach_state_', num2str(BFields(i)), 'G_', basis, '.mat']);
    cAll{i} = load(['data/c3Sigma_state_', num2str(BFields(i)), 'G_', basis, '.mat']);
end


power_up = 1; % W
power_dn = 1; % W

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


for i = 1:length(BFields)
    %% ignore quantum numbers that are not common to all 3 states
    f = fAll{i};
    c = cAll{i};
    common_qnums = intersect(c.qnums.Properties.VariableNames,f.qnums.Properties.VariableNames);

    f.qnums(:,~ismember(f.qnums.Properties.VariableNames,common_qnums)) = [];
    c.qnums(:,~ismember(c.qnums.Properties.VariableNames,common_qnums)) = [];

    % ignore some quantum numbers in rotational matrix element
    qnums_ignore = {'S','Sigma','Lambda'};
    f_ignore = ismember(f.qnums.Properties.VariableNames,qnums_ignore);
    c_ignore = ismember(c.qnums.Properties.VariableNames,qnums_ignore);

    %% find energy reference point
    switch basis
        case 'aFC'
            [~,ind] = evec_ind({'J','I','F','m_F'},[1,5,6,4],c,c.psi);
        case 'aUC'
            [~,ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],c,c.psi);
        case 'aIC'
            [~,ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
    end
    E0_up = mean(c.E(ind)-min(f.E));

    %% transition dipole moments, angular momentum part
    p = -1:1; % spherical index

    T1q_c_f = sphten([0 0 1]);

    % transition dipole moments
    switch basis
        case 'aFC'
            rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_aFC,...
                {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_c_f);
        case {'aUC','aIC'}
            rot_TDM_c_f_components = operator_matrix(@transition_dipole_case_a,...
                {c.qnums(:,~c_ignore),f.qnums(:,~f_ignore)},{'eta','J','Omega','m_J'},p,T1q_c_f);
    end

    % off-diagonal Hamiltonian blocks
    rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T_up(p+2),1,1,3),3);

    %% transition dipole moments, vibronic part
    a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
    X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

    a3S_c3S_elecTDM = (sqrt(1-const.c3Sigma.singlet_fraction)*(f.qnums.S==1)...
        .*interp1(a3S_c3S_elecTDM_data(:,1)*const.abohr,a3S_c3S_elecTDM_data(:,2)*const.e*const.abohr,c.r) ...
        + sqrt(const.c3Sigma.singlet_fraction)*(f.qnums.S==0)...
        .*interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,c.r));
    a3S_c3S_elecTDM(isnan(a3S_c3S_elecTDM)) = 0;
    rovib_TDM_c_f = rot_TDM_c_f.*permute(a3S_c3S_elecTDM,[3 1 2]);

%     X1S_c3S_elecTDM_interp = sqrt(const.c3Sigma.singlet_fraction)...
%         *interp1(X1S_B1P_elecTDM_data(:,1)*const.abohr,X1S_B1P_elecTDM_data(:,2)*const.e*const.abohr,X.r,'spline');

    %% laser hamiltonian matrix elements
    psi_c = c.psi.*permute(c.psi_r,[1 3 2]);
    psi_feshbach_interp = permute(interp1(f.r,f.psi(:,:,1)',c.r,'spline')',[1 3 2]);
    TDM_R_c_f = mtimesx(conj(permute(psi_c,[2 1 3])),mtimesx(rovib_TDM_c_f,psi_feshbach_interp));
    TDM_c_f = trapz(c.r,TDM_R_c_f,3);
    H_up = TDM_c_f*Efield_up;

%     psi_c_interp = interp1(c.r,c.psi_r,X.r,'spline');
%     vibronic_TDM_X_c = trapz(X.r,X.psi_r.*X1S_c3S_elecTDM_interp.*psi_c_interp);
%     H_dn = vibronic_TDM_X_c * Efield_dn * (c.psi'*rot_TDM_X_c*X.psi);

    % rabi frequencies per sqrt(mW)
    switch basis
        case 'aFC'
            [q_ind,v_ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],c,c.psi);
            mj0Ind = find(12 < c.E*1e-9/const.h - 325100 & c.E*1e-9/const.h - 325100 < 13);
        case 'aUC'
            [q_ind,v_ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,0,3/2,5/2],c,c.psi);
            mj0Ind = find(12 < c.E*1e-9/const.h - 325100 & c.E*1e-9/const.h - 325100 < 13);
        case 'aIC'
            [q_ind,v_ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
    end
    tdmUp(i) = max(abs(H_up(v_ind)));
    rabi_up_sqrt_mW = tdmUp(i)/const.h*1e-6 *sqrt(1e-3/power_up);
    fprintf('up leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_up_sqrt_mW)

    tdmUpSum(i) = sqrt(sum(H_up(mj0Ind).^2));
    
%     switch basis
%         case 'aFC'
%             [~,v_ind2] = evec_ind({'J','I','F','m_F'},[0,5,5,4],X,X.psi);
%         case 'aIC'
%             [~,v_ind2] = evec_ind({'J','m_J','I','m_I'},[0,0,5,4],X,X.psi);
%         case 'aUC'
%             [~,v_ind2] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[0,0,3/2,5/2],X,X.psi);
%     end
%     rabi_dn_sqrt_mW = max(abs(H_dn(v_ind,v_ind2)))/const.h*1e-6 *sqrt(1e-3/power_dn);
%     fprintf('down leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_dn_sqrt_mW)

end


