
% Calculate rotational and hyperfine structure of c3Sigma state of NaCs

clear;

c = constants();

Jmax = 3;
mtot = [3 5];
B = 855*1e-4;

%% build basis
basis.aUC.qnums = build_basis({'Lambda','Omega','J','S','i_Na','i_Cs'},...
    {c.c3Sigma.Lambda,c.c3Sigma.Omega,c.c3Sigma.Omega:Jmax,c.c3Sigma.S,c.i_Na,c.i_Cs},[2 2 2 2 1 1],'a');
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

basis.aUC.ops.Fz = basis.aUC.ops.J_z+basis.aUC.ops.i_Na_z+basis.aUC.ops.i_Cs_z;

basis.aUC.ops.H = basis.aUC.ops.Hrot + basis.aUC.ops.Hsrot ...
    + basis.aUC.ops.HZ0elecspin.*B + basis.aUC.ops.Hhf_Na + basis.aUC.ops.Hhf_Cs + basis.aUC.ops.H_OmegaDoubling;

%% truncate basis
[basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.J_z+ops.i_Na_z+ops.i_Cs_z,mtot);

%% diagonalize
[evecs,evals] = eig(basis.aUC.ops.H);
evals = real(diag(evals));

%% plot
mF_ev = real(diag(evecs'*basis.aUC.ops.Fz*evecs));
figure(1);
clf;
hold on; box on;
for i = 1:numel(evals)
    plot(mF_ev(i) + [-0.4 0.4],[1 1]*evals(i)*1e-9/c.h,'-k')
end
hold off;
xlabel('M_{tot}');
ylabel('Energy (GHz)');

%% solve the vibrational problem
% simulation parameters
Nx = 2000;
rmin = 5;
rmax = 100;

Erange = 0.049415 + [-1 1]*1e-4; 

% total hamiltonian
W =@(r) NaCscPES(r);

% call the solver
[E_out,nodes_out,err_est,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,c.mu_nacs/c.me,1,1);

figure(2);
clf;
for i = 1:size(psi,3)
    subplot(size(psi,3),1,i)
    plot(r,psi(:,:,i),'linewidth',2);
    set(gca,'xscale','log')
    xlim([min(r) max(r)])
    ylabel('\psi(R)')
end
xlabel('R (a_0)')


%% save data
qnums = basis.aUC.qnums;
psi = evecs;
E = evals;
save(['c3Sigma_state_' num2str(round(B*1e4)) 'G.mat'],'qnums','psi','E');


