clear;
c = constants();

% Coupled-channel calculation in ground state of NaCs. Configured to
% find feshbach state and lowest 2 trap states.

% simulation parameters
Nx = 8000;
rmin = 4.5;
rmax = 10000;
mtot = 4; % total angular momentum to truncate basis
B = (855)*1e-4; %/c.B_au; % number in parens is in gauss

% uncoupled basis
basis.UC = join_basis(atom_basis('Na'),atom_basis('Cs'));
basis_mol.qnums = build_basis({'eta','Lambda','N'},{1,0,0},[0,60,1],'b');
basis_mol.ops = build_operators(basis_mol.qnums);
basis.UC = join_basis(basis.UC,basis_mol);

% spin-coupled basis
[basis.SC,basis.UC,basis_change.UC_SC] = couple_angmom(basis.UC,'s_Na','s_Cs','S');
basis = rmfield(basis,'UC');

% projectors onto singlet and triplet subspaces
basis.SC.ops.Proj_singlet = 1*diag(basis.SC.qnums.S==0);
basis.SC.ops.Proj_triplet = 1*diag(basis.SC.qnums.S==1);

% atomic hamiltonians
basis.SC.ops.H0 = (basis.SC.ops.H_Na_hyperfine + basis.SC.ops.H_Cs_hyperfine ...
    + (basis.SC.ops.H0_Na_zeeman + basis.SC.ops.H0_Cs_zeeman)*B)/c.hartree;

% form hund's case b basis
[basis.b,basis.SC,basis_change.SC_b] = couple_angmom(basis.SC,'N','S','J');

% transformation to hund's case a
basis.a.qnums = build_basis({'eta','s_Na','s_Cs','i_Na','i_Cs','J','S','Lambda'},...
    {1,c.s_Na,c.s_Cs,c.i_Na,c.i_Cs,0:1,0:1,0},[0 0 0 1 1 0 0 0],'a');
basis.a.qnums((basis.a.qnums.J==1 & basis.a.qnums.S==0) | (basis.a.qnums.J==0 & basis.a.qnums.S==1),:) = [];
basis.a.ops = build_operators(basis.a.qnums);

% trap
waist = 1064e-9; % meter
Power = 10e-3; % watt
I0 = 2*Power/(pi*waist^2); % W/m^2
trap_depth = (c.apol_Cs*I0)/(2*c.eps0*c.c)/c.hartree;
atrap = 2*((c.abohr)/waist).^2;
Vtrap =@(r) trap_depth*(1-exp(-atrap*r.^2));

% spin dipole-dipole interactions give rise to the FB resonance
% include electron-electron spin and electron-nuclear spin
basis.SC.ops.Vdd0 = (c.mu0*c.uB^2/(4*pi*c.hartree*c.abohr^3))*(...
    c.gs_Na*c.gs_Cs*dd_operator(basis.SC.ops,'s_Na','s_Cs') ...
    + c.gs_Na*c.gi_Cs*dd_operator(basis.SC.ops,'s_Na','i_Cs') ...
    + c.gi_Na*c.gs_Cs*dd_operator(basis.SC.ops,'i_Na','s_Cs'));

% truncate basis to specific mtot
[basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.J_z+ops.i_Na_z+ops.i_Cs_z,mtot);

% total hamiltonian
W =@(r) herm(basis.SC.ops.H0 + NaCsXPES(r).*basis.SC.ops.Proj_singlet ...
    + NaCsaPES(r).*basis.SC.ops.Proj_triplet + basis.SC.ops.Vdd0./(r.^3) ...
    + Vtrap(r).*basis.SC.ops.I);


% find lowest energy channel and set energy search range
rtest = logspace(log10(rmin),log10(rmax),1e2);
Wtest = W(reshape(rtest,1,1,[]));
[V,D] = eigenshuffle(Wtest);
r_threshold = 300;
[E_lowest_chan_threshold,ind] = min(D(:,abs(rtest-r_threshold)==min(abs(rtest-r_threshold))));
Erange = E_lowest_chan_threshold + [-50e6 1.2e6]*c.h/c.hartree;

figure(1);
clf;
plot(rtest,D'*c.hartree/c.h*1e-9,'linewidth',2);
ylim([-1 1]*20)
set(gca,'xscale','log','fontsize',14);
xlabel('R (a_0)');
ylabel('E (GHz)');
drawnow();

% call the solver
[E_out,nodes_out,err_est,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,c.mu_nacs/c.me);

% display results
diff(E_out)*c.hartree/c.h * 1e-6
errstr(E_out/c.wavenum2hartree*(29.9792458),err_est/c.wavenum2hartree*(29.9792458))

figure(2);
clf;
for i = 1:size(psi,3)
    subplot(size(psi,3),1,i)
    plot(r,psi(:,:,i),'linewidth',2),set(gca,'xscale','log');
    set(gca,'fontsize',18)
    xlim([min(r) max(r)])
    ylabel('\psi(R)')
end
xlabel('R (a_0)')

% change basis to hund's case (a) for calculating electric dipole
% transition matrix elements
b_a = operator_matrix(@case_b2a_element,{basis.a.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});
psi_a = zeros(size(basis.a.qnums,1),numel(r),size(basis.b.qnums,1));
for i = 1:size(psi,3)
    psi_a(:,:,i) = b_a*psi(:,:,i);
end

% save data
qnums = basis.a.qnums;
r = r*c.abohr;
psi = psi_a;
E = E_out*c.hartree;
save(['feshbach_state_' num2str(round(B*1e4)) 'G.mat'],'qnums','r','psi','E','B');
