%Calculate upleg effective Hamiltonian for different c3Sigma parameters
%with most of the computational heavy lifting already done

function pre = precomp_ops(r)

pre.basis.aUC.qnums = build_basis({'eta','Lambda','Omega','J','S','i_Na','i_Cs'},...
    {2,r.c3Sigma.Lambda,r.c3Sigma.Omega,r.c3Sigma.Omega:r.Jmax,r.c3Sigma.S,r.i_Na,r.i_Cs},[0 2 2 2 2 1 1],'a');
pre.basis.aUC.ops = build_operators(pre.basis.aUC.qnums);

Sp = operator_matrix(@Sp_case_a_element,pre.basis.aUC.qnums,{'J','Omega','m_J','S'},1);
Sm = operator_matrix(@Sp_case_a_element,pre.basis.aUC.qnums,{'J','Omega','m_J','S'},-1);
pre.basis.aUC.ops.S_x = (Sm-Sp)/sqrt(2);
pre.basis.aUC.ops.S_y = -(Sm+Sp)/(sqrt(2)*1i);
pre.basis.aUC.ops.S_z = operator_matrix(@Sp_case_a_element,pre.basis.aUC.qnums,{'J','Omega','m_J','S'},0);

pre.basis.aUC.ops.omdub = operator_matrix(@c3Sigma_omegadoubling_element,pre.basis.aUC.qnums,{'J','Omega','Lambda','Sigma'});
pre.basis.aUC.ops.szel = operator_matrix(@Sz_case_a_element,pre.basis.aUC.qnums,{'J','Omega','m_J','S'});
pre.basis.aUC.ops.iNa = op_dot(pre.basis.aUC.ops,'i_Na','S');
pre.basis.aUC.ops.iCs = op_dot(pre.basis.aUC.ops,'i_Cs','S');

[pre.basis.aIC,pre.basis.aUC,pre.basis.change.aUC_IC] = couple_angmom(pre.basis.aUC,'i_Na','i_Cs','I');
[pre.basis.aFC,pre.basis.aIC,pre.basis.change.IC_FC] = couple_angmom(pre.basis.aIC,'J','I','F');
pre.basis.change.aUC_FC = pre.basis.change.aUC_IC*pre.basis.change.IC_FC;

pre.c3Sigma = c3Sigma(r.B,'aFC',0,r.Jmax,r.mtot,r.f_vib,r);
if r.N_tot > 0
    pre.f = feshbach_rot(r.B,'aFC',1,r,r.N_tot);
else
    pre.f = feshbach(r.B,'aFC',1,r);
end
pre.common_qnums = intersect(pre.c3Sigma.qnums.Properties.VariableNames,pre.f.qnums.Properties.VariableNames);

% ignore qnums not common to all bases
pre.f.qnums(:,~ismember(pre.f.qnums.Properties.VariableNames,pre.common_qnums)) = [];
pre.c3Sigma.qnums(:,~ismember(pre.c3Sigma.qnums.Properties.VariableNames,pre.common_qnums)) = [];

% ignore hund's case (a) quantum numbers in rotational matrix element --
% use only hund's case (c)
qnums_ignore = {'S','Sigma','Lambda'};
pre.f_ignore = ismember(pre.f.qnums.Properties.VariableNames,qnums_ignore);
pre.c_ignore = ismember(pre.c3Sigma.qnums.Properties.VariableNames,qnums_ignore);

p = -1:1; % spherical index
T1q_f_z = sphten([0 0 1]);
T1q_f_perp = sphten([1 1 0]/sqrt(2));

rot_TDM_f_z_components = operator_matrix(@transition_dipole_case_aFC,...
    {pre.c3Sigma.qnums(:,~pre.c_ignore),pre.f.qnums(:,~pre.f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_f_z);
rot_TDM_f_perp_components = operator_matrix(@transition_dipole_case_aFC,...
    {pre.c3Sigma.qnums(:,~pre.c_ignore),pre.f.qnums(:,~pre.f_ignore)},{'eta','J','Omega','I','F','m_F'},p,T1q_f_perp);

rot_TDM_f_z = (-1).^reshape(p,1,1,3) .* rot_TDM_f_z_components;
rot_TDM_f_perp = (-1).^reshape(p,1,1,3) .* rot_TDM_f_perp_components;

% push polarization info to 4th dimension
pre.rot_TDM_f_z = permute(rot_TDM_f_z,[1,2,4,3,5]);
pre.rot_TDM_f_perp = permute(rot_TDM_f_perp,[1,2,4,3,5]);

a3S_b3P_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_b3P_Full');
a3S_c3S_elecTDM_data = importdata('lib/rosario_potentials/DM_a3S_c3S_Full');
X1S_B1P_elecTDM_data = importdata('lib/rosario_potentials/DM_X1S_B1P_Full');

X1S_B1P_elecTDM_interp2 = interp1(X1S_B1P_elecTDM_data(:,1)*r.abohr,X1S_B1P_elecTDM_data(:,2)*r.e*r.abohr,pre.c3Sigma.r,'spline');
a3S_c3S_elecTDM_interp = interp1(a3S_c3S_elecTDM_data(:,1)*r.abohr,a3S_c3S_elecTDM_data(:,2)*r.e*r.abohr,pre.c3Sigma.r,'spline');
a3S_b3P_elecTDM_interp = interp1(a3S_b3P_elecTDM_data(:,1)*r.abohr,a3S_b3P_elecTDM_data(:,2)*r.e*r.abohr,pre.c3Sigma.r,'spline');

%% laser hamiltonian matrix elements

% interpolate radial coord
psi_feshbach_interp = permute(interp1(pre.f.r,permute(pre.f.psi,[2 1 3]),pre.c3Sigma.r,'spline'),[2 1 3]);

pre.f_X_rows = all(pre.f.qnums{:,{'Lambda','S'}}==[0 0],2);
pre.f_a_rows = all(pre.f.qnums{:,{'Lambda','S'}}==[0 1],2);
pre.c_B_rows = strcmp(cellstr(pre.c3Sigma.qnums_vib.term),'B1P1');
pre.c_b_rows = strcmp(cellstr(pre.c3Sigma.qnums_vib.term),'b3P1');
pre.c_c_rows = strcmp(cellstr(pre.c3Sigma.qnums_vib.term),'c3S1');

pre.vibronic_TDM_B_f = trapz(pre.c3Sigma.r,psi_feshbach_interp(pre.f_X_rows,:,:).*X1S_B1P_elecTDM_interp2.*pre.c3Sigma.psi_vib(pre.c_B_rows,:),2);
pre.vibronic_TDM_b_f = trapz(pre.c3Sigma.r,psi_feshbach_interp(pre.f_a_rows,:,:).*a3S_b3P_elecTDM_interp.*pre.c3Sigma.psi_vib(pre.c_b_rows,:),2);
pre.vibronic_TDM_c_f = trapz(pre.c3Sigma.r,psi_feshbach_interp(pre.f_a_rows,:,:).*a3S_c3S_elecTDM_interp.*pre.c3Sigma.psi_vib(pre.c_c_rows,:),2);

pre.psi_init = zeros(size(pre.f.E));
if r.initState == 1
    indState = find(pre.f.E < pre.f.E_lowest_chan_threshold);
    pre.psi_init(indState(end)) = 1; %Find highest energy state below threshold
elseif r.initState == 2
    indState = find(pre.f.E > pre.f.E_lowest_chan_threshold);
    pre.psi_init(indState(1)) = 1; %Find lowest energy state above threshold
end

end