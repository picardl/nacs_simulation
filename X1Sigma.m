
clear;


const = constants();

Nmax = 0;
gamma_pol = 0;

% simulation parameters
Nx = 500;
rmin = 4;
rmax = 15;
Erange = -0.0224 + [-1 1]*1e-4;

B = 855e-4;

%% build operators in uncoupled basis
bases = {'UC','IC','FC','F1C','F2C'};

basis.UC.qnums = build_basis({'eta','s_Na','s_Cs','Lambda','i_Na','i_Cs','N','S'},...
    {3,const.s_Na,const.s_Cs,0,const.i_Na,const.i_Cs,0:Nmax,0},[0 0 0 0 1 1 1 1]);
basis.UC.ops = build_operators(basis.UC.qnums);
Nstates = size(basis.UC.qnums,1);

%% successively couple angular momenta to form coupled bases
[basis.IC,basis.UC,basis.change.UC_IC] = couple_angmom(basis.UC,'i_Na','i_Cs','I');
[basis.FC,basis.IC,basis.change.IC_FC] = couple_angmom(basis.IC,'I','N','F');
basis.change.UC_FC = basis.change.UC_IC*basis.change.IC_FC;
[basis.F1C,~,basis.change.UC_F1C] = couple_angmom(basis.UC,'i_Na','N','F1');
[basis.F2C,~,basis.change.UC_F2C] = couple_angmom(basis.UC,'i_Cs','N','F2');
[basis.b,basis.UC,basis.change.UC_b] = couple_angmom(basis.UC,'N','S','J');

%% spin-spin scalar coupling, i_Na.i_Cs
% do it in all bases to check consistency
for a = bases
    basis.(a{1}).ops.i_Nadoti_Cs = op_dot(basis.(a{1}).ops,'i_Na','i_Cs');
end

%% spin-rotation scalar coupling, i_Na.N and i_Cs.N
basis.FC.ops.i_NadotN = op_dot(basis.FC.ops,'i_Na','i_Cs');
for a = bases
    basis.(a{1}).ops.i_NadotN = op_dot(basis.(a{1}).ops,'i_Na','N');
    basis.(a{1}).ops.i_CsdotN = op_dot(basis.(a{1}).ops,'i_Cs','N');
end

%% electric quadrupole coupling
j1_dot_j2_cpl =@(f,j1,j2) (f.*(f+1)-j1.*(j1+1)-j2.*(j2+1))/2;
EQ =@(F,N,I) -(3*j1_dot_j2_cpl(F,N,I).^2 + 3/2*j1_dot_j2_cpl(F,N,I) - I.*(I+1).*N.*(N+1))./(2*I.*(2*I-1).*(2*N-1).*(2*N+3));
basis.F1C.ops.EQ1 = diag(EQ(basis.F1C.qnums.F1,basis.F1C.qnums.N,basis.F1C.qnums.i_Na));
basis.F2C.ops.EQ2 = diag(EQ(basis.F2C.qnums.F2,basis.F2C.qnums.N,basis.F2C.qnums.i_Cs));

%% convert all operators back to uncoupled basis for all subsequent calculations
UC_ops = fields(basis.UC.ops);
for a = 2:numel(bases)
    b_ops = fields(basis.(bases{a}).ops);
    not_in_UC = find(~ismember(b_ops,UC_ops))';
    for i = not_in_UC
        basis.UC.ops.(b_ops{i}) = basis.change.(['UC_' bases{a}])*basis.(bases{a}).ops.(b_ops{i})*basis.change.(['UC_' bases{a}])';
    end
end
basis = rmfield(basis,{'IC','FC','F1C','F2C'});


%% transformation to hund's case a -- turns out it's the identity
basis.a.qnums = build_basis({'eta','i_Na','i_Cs','J','S','Lambda'},...
    {3,const.i_Na,const.i_Cs,0:Nmax,0,0},[0 1 1 2 0 0],'a');
basis.a.ops = struct();
basis.change.b_a = operator_matrix(@case_b2a_element,{basis.a.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});
basis.change.UC_a = basis.change.b_a*basis.change.UC_b;
f = fields(basis.UC.ops);
for i = 1:numel(f)
    basis.a.ops.(f{i}) = basis.change.UC_a*basis.UC.ops.(f{i})*basis.change.UC_a';
end

%% truncate basis
[basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.J_z+ops.i_Na_z+ops.i_Cs_z,[2 3 4 5]);

%% build hamiltonian
basis.UC.ops.H0 = const.X1Sigma.Bv * basis.UC.ops.N_sq + ...
    const.X1Sigma.c4 * basis.UC.ops.i_Nadoti_Cs + ...
    const.X1Sigma.c1 * basis.UC.ops.i_NadotN + ...
    const.X1Sigma.c2 * basis.UC.ops.i_CsdotN + ...
    const.X1Sigma.eQq1 * basis.UC.ops.EQ1 + ...
    const.X1Sigma.eQq2 * basis.UC.ops.EQ2;
basis.UC.ops.Hz0 = -const.X1Sigma.gr*const.uN*( basis.UC.ops.N_z ) + ...
    -const.X1Sigma.g1*const.uN*(1-const.X1Sigma.sigma1)*( basis.UC.ops.i_Na_z ) + ...
    -const.X1Sigma.g2*const.uN*(1-const.X1Sigma.sigma2)*( basis.UC.ops.i_Cs_z );

H = basis.UC.ops.H0 + basis.UC.ops.Hz0.*B;

[psi,E] = eigenshuffle(H);
E = real(E);

mF = diag(psi'*basis.UC.ops.F_z*psi);
figure(1);
clf;
hold on; box on;
for i = 1:numel(E)
    plot([-1 1]*0.4+mF(i),[1 1]*E(i)/const.h*1e-6,'-k');
end
hold off;
xlabel('m_F');
ylabel('E (MHz)');

%% solve the vibrational problem
rtest = linspace(rmin,rmax,1e3);

% total hamiltonian
W =@(r) NaCsXPES(r);
Wtest = W(rtest);

% call the solver
[E_vib,nodes_out,err_est,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,@NaCsXPES,Erange,const.mu_nacs/const.me,1,1);

figure(2);
clf;
hold on;
box on;
for i = 1:size(psi_r,3)
    subplot(size(psi_r,3),1,i)
    plot(rtest,Wtest);
    plot([rmin rmax],min(Erange)*[1 1],'-k');
    plot([rmin rmax],max(Erange)*[1 1],'-k');
    plot(r,psi_r(:,:,i)*0.1*peak2peak(Wtest) + E_vib,'linewidth',2);
    set(gca,'xscale','log')
    xlim([min(r) max(r)])
    ylabel('\psi(R)')
end
xlabel('R (a_0)')
hold off;

%% save data
psi = basis.change.UC_a*psi;
qnums = basis.a.qnums;
E = E + E_vib*const.hartree;

save(['X1Sigma_state_' num2str(round(B*1e4)) 'G.mat'],'qnums','psi','psi_r','r','E','B');
