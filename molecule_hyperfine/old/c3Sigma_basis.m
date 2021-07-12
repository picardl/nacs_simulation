
% Calculate rotational and hyperfine structure of c3Sigma state of NaCs

clear;

c = constants();

Jmax = 3;

save_basis = 'aUC';

%% build bases
basis.aUC.qnums = build_basis({'eta','Lambda','Omega','J','S','i_Na','i_Cs'},...
    {2,c.c3Sigma.Lambda,c.c3Sigma.Omega,c.c3Sigma.Omega:Jmax,c.c3Sigma.S,c.i_Na,c.i_Cs},[0 2 2 2 2 1 1],'a');
basis.aUC.ops = build_operators(basis.aUC.qnums);

%% build hamiltonian
basis.aUC.ops.Hrot = basis.aUC.ops.J_sq;

basis.aUC.ops.H_OmegaDoubling = operator_matrix(@c3Sigma_omegadoubling_element,basis.aUC.qnums,{'J','Omega','Lambda','Sigma'});

basis.aUC.ops.Hsrot = operator_matrix(@spin_rotation_element,basis.aUC.qnums,{'J','Omega','m_J','S'});

basis.aUC.ops.HZ0elecspin = c.uB*operator_matrix(@Sz_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'});

Sp = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},1);
Sm = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},-1);
basis.aUC.ops.S_x = (Sm-Sp)/sqrt(2);
basis.aUC.ops.S_y = -(Sm+Sp)/(sqrt(2)*1i);
basis.aUC.ops.S_z = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},0);

basis.aUC.ops.Hhf_Na = op_dot(basis.aUC.ops,'i_Na','S');
basis.aUC.ops.Hhf_Cs = op_dot(basis.aUC.ops,'i_Cs','S');

%% go to fully coupled basis
[basis.aIC,basis.aUC,basis.change.aUC_IC] = couple_angmom(basis.aUC,'i_Na','i_Cs','I');
[basis.aFC,basis.aIC,basis.change.IC_FC] = couple_angmom(basis.aIC,'J','I','F');
basis.change.aUC_FC = basis.change.aUC_IC*basis.change.IC_FC;

%% transform F operators back to uncoupled basis
f = setdiff(fields(basis.aFC.ops),fields(basis.aUC.ops));
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = basis.change.aUC_FC*basis.aFC.ops.(f{i})*basis.change.aUC_FC';
end

[basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.F_z,[2 3 4 5]);

%% solve the vibrational problem
% simulation parameters
Nx = 1000;
rmin = 5;
rmax = 20;

Erange = 0.049415 + [-1 1]*1e-4; 

% call the solver
[E_vib,nodes_out,err_est,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,@NaCscPES,Erange,c.mu_nacs/c.me,1,1);

%% save
qnums = basis.(save_basis).qnums;
ops = basis.(save_basis).ops;
r = r*c.abohr;
psi_r = psi_r/sqrt(c.abohr);

fname = '../data/c3Sigma_basis.mat';
save(fname,'basis','E_vib','r','psi_r','qnums','ops');
disp(['saved ' fname]);