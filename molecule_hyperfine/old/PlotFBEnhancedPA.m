clear;

const = constants();

% Determine tdm dependence on B field close to FB resonance
% Needs to be in nacs_simulation/molecule_hyperfine folder! 

basis = 'aUC';
%% load data
BFields = [850:1:860];
tdmUp = zeros(size(BFields));
tdmUpSum = zeros(size(BFields));
tdmUpMax = zeros(size(BFields));
rabi_up_sqrt_mW = zeros(size(BFields));

fAll = {};
cAll = {};
for i = 1:length(BFields)
    fAll{i} = load(['data/fb_', num2str(BFields(i)), 'G_', basis, '_220124.mat']);
    cAll{i} = load(['../data/c_' [strrep(num2str(BFields(i)),'.','p') 'G'] '_vib325THz_' basis '_220124.mat']);
end


power_up = 1; % W
power_dn = 1; % W

waist = 13e-6; % m

pol_up = 'sigp';
pol_dn = 'sigp';

%% laser stuff
% electric fields
Efield_up = sqrt(4*const.eta0*power_up/(pi*waist^2)); % V/m
Efield_dn = sqrt(4*const.eta0*power_dn/(pi*waist^2)); % V/m

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

%%
%Find and plot TDMs

for i = 1:length(BFields)
    %% ignore quantum numbers that are not common to all 3 states
    f = fAll{i}.out;
    c = cAll{i}.out;
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
    psi_c = c.psi.*permute(c.r,[1 3 2]);
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
            [q_ind,v_ind] = evec_ind({'J','I','F','m_F'},[1,5,6,4],c,c.psi);
            mj0Ind = find(12 < c.E*1e-9/const.h - 325100 & c.E*1e-9/const.h - 325100 < 13);
        case 'aUC'
            [q_ind,v_ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,0,3/2,5/2],c,c.psi);
            mj0Ind = find(12 < c.E*1e-9/const.h - 325100 & c.E*1e-9/const.h - 325100 < 13);
        case 'aIC'
            [q_ind,v_ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
    end
    tdmUp(i) = max(abs(H_up(v_ind)));


    tdmUpSum(i) = sqrt(sum(H_up(mj0Ind).^2));
    tdmUpMax(i) = max(max(abs(H_up(mj0Ind))));
    rabi_up_sqrt_mW(i) = tdmUpMax(i)*sqrt(1e-3/power_up);
    fprintf('up leg transition strength = %1.3g Hz/sqrt(mW)\n',rabi_up_sqrt_mW)
    
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

%%
%Plotting

sg = load('C:\NaCs1pt5_Data_temp\Data\20200724\data_20200724_182304.mat');
surv = sg.Analysis.SurvivalProbability(3,:);

Gamma = 100e6*const.h;
pow = 3.7;
rabi_up = rabi_up_sqrt_mW*sqrt(pow);
t = 2e-4;
Rscatter_atom = (Gamma/(2*const.hbar)) * 2*abs(rabi_up).^2./(2*abs(rabi_up).^2 + Gamma.^2);
popn = 0.8*exp(-Rscatter_atom*t);

fitResult = fminsearch(@(x) sum((x(2)*exp(-Rscatter_atom*x(1)) - surv).^2),[2e-4,0.8]);


% dataObj = NaCsData(20200724,182304);
% [allChnsSurvival, allChnsSurvivalErr,allScanParams,plotted] = plotSurvivalAvg(dataObj,0,{[1,2]}, {[3,4]},0,0);
% set(plotted{1},'LineStyle','None')
% set(plotted{1},'MarkerSize',6)
% set(plotted{1},'MarkerFaceColor',[0.00,0.58,0.96])
% set(plotted{1},'XData',[-3:0.5:3])
% ax = gca();
% set(ax,'Title',[])
% ylabel('Na + Cs Survival')
% xlabel('B - B_{FR} [G]', 'Interpreter','tex')
% grid off

% hold on
popn = fitResult(2)*exp(-Rscatter_atom*fitResult(1));
simDat = plot(BFields + 10.5 - 866, popn);
set(simDat,'LineStyle','--')
set(simDat,'Color',[0,0,0])
set(simDat,'LineWidth',2)

legend(['Experimental Data','Effective Hamiltonian Model'])



