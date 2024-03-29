
% Calculate rotational and hyperfine structure of c3Sigma state of NaCs

clear basis;

c = constants();

Jmax = 3;
mtot = [2 3 4 5]; % [3 4 5];
% B = 853*1e-4;
save_basis = 'aUC';

%% build bases
basis.aUC.qnums = build_basis({'eta','Lambda','Omega','J','S','i_Na','i_Cs'},...
    {2,c.c3Sigma.Lambda,c.c3Sigma.Omega,c.c3Sigma.Omega:Jmax,c.c3Sigma.S,c.i_Na,c.i_Cs},[0 2 2 2 2 1 1],'a');
basis.aUC.ops = build_operators(basis.aUC.qnums);

%% build hamiltonian
basis.aUC.ops.Hrot = c.c3Sigma.Be*basis.aUC.ops.J_sq;

basis.aUC.ops.H_OmegaDoubling = c.c3Sigma.wef*operator_matrix(@c3Sigma_omegadoubling_element,basis.aUC.qnums,{'J','Omega','Lambda','Sigma'});

basis.aUC.ops.Hsrot = (c.c3Sigma.gamma-2*c.c3Sigma.Be).*operator_matrix(@spin_rotation_element,basis.aUC.qnums,{'J','Omega','m_J','S'});

basis.aUC.ops.HZ0elecspin = c.c3Sigma.gS*c.uB*operator_matrix(@Sz_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'});

Sp = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},1);
Sm = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},-1);
basis.aUC.ops.S_x = (Sm-Sp)/sqrt(2);
basis.aUC.ops.S_y = -(Sm+Sp)/(sqrt(2)*1i);
basis.aUC.ops.S_z = operator_matrix(@Sp_case_a_element,basis.aUC.qnums,{'J','Omega','m_J','S'},0);

basis.aUC.ops.Hhf_Na = c.c3Sigma.alpha1*op_dot(basis.aUC.ops,'i_Na','S');
basis.aUC.ops.Hhf_Cs = c.c3Sigma.alpha2*op_dot(basis.aUC.ops,'i_Cs','S');

basis.aUC.ops.H = basis.aUC.ops.Hrot + basis.aUC.ops.Hsrot ...
    + basis.aUC.ops.HZ0elecspin.*B + basis.aUC.ops.Hhf_Na + basis.aUC.ops.Hhf_Cs + basis.aUC.ops.H_OmegaDoubling;

%% go to fully coupled basis
[basis.aIC,basis.aUC,basis.change.aUC_IC] = couple_angmom(basis.aUC,'i_Na','i_Cs','I');
[basis.aFC,basis.aIC,basis.change.IC_FC] = couple_angmom(basis.aIC,'J','I','F');
basis.change.aUC_FC = basis.change.aUC_IC*basis.change.IC_FC;

%% transform F operators back to uncoupled basis
f = setdiff(fields(basis.aFC.ops),fields(basis.aUC.ops));
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = basis.change.aUC_FC*basis.aFC.ops.(f{i})*basis.change.aUC_FC';
end
% max(reshape(abs(basis.aUC.ops.J_z + basis.aUC.ops.i_Na_z + basis.aUC.ops.i_Cs_z - basis.aUC.ops.F_z),[],1))

%% truncate basis
% basis = rmfield(basis,'aUC');
[basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.F_z,mtot);

%% diagonalize
[evecs,evals] = eig(basis.(save_basis).ops.H);
evals = abs(diag(evals));

%% plot
% mF_expect = real(diag(evecs'*basis.(save_basis).ops.F_z*evecs));
% figure(1);
% clf;
% hold on; box on;
% for i = 1:numel(evals)
%     plot(mF_expect(i) + [-0.4 0.4],[1 1]*evals(i)*1e-9/c.h,'-k')
% end
% hold off;
% xlabel('M_{tot}');
% ylabel('Energy (GHz)');

%% solve the vibrational problem
% simulation parameters
Nx = 1000;
rmin = 5;
rmax = 20;

Erange = 0.049415 + [-1 1]*1e-4; 

% call the solver
[E_vib,nodes_out,err_est,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,@NaCscPES,Erange,c.mu_nacs/c.me,1,1);

% figure(2);
% clf;
% for i = 1:size(psi_r,3)
%     subplot(size(psi_r,3),1,i)
%     plot(r,psi_r(:,:,i),'linewidth',2);
%     set(gca,'xscale','log')
%     xlim([min(r) max(r)])
%     ylabel('\psi(R)')
% end
% xlabel('R (a_0)')

%% save data
qnums = basis.(save_basis).qnums;
ops = basis.(save_basis).ops;
psi = evecs;
r = r*c.abohr;
psi_r = psi_r/sqrt(c.abohr);
E = evals + E_vib*c.hartree;
fname = ['data/c3Sigma_state_' num2str(B*1e4) 'G_' save_basis '.mat'];
save(fname,'qnums','ops','psi','psi_r','r','E','B');
disp(['saved file ' fname])
