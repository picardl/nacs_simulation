function out = Xa_complex

c = constants();

save_basis = 'aUC';

% Coupled-channel calculation in ground state of NaCs. 

%% simulation parameters
Nx = 2000;
rmin = 4.5;
rmax = 1e2;
mtot = 4; % total angular momentum to truncate basis
Erange = [-0.023 -1e-4];
% Erange_vs_threshold = [-0.023 0];

%% uncoupled basis
basis.UC.qnums = build_basis({'s_Na','s_Cs'},{1/2,1/2},[1,1]);
basis.UC.ops = build_operators(basis.UC.qnums);
basis_mol.qnums = build_basis({'eta','Lambda','N'},{1,0,0:2},[0,0,1],'b');
basis_mol.ops = build_operators(basis_mol.qnums);
basis.UC = join_basis(basis.UC,basis_mol);

[basis.SC,basis.UC,basis.change.UC_SC] = couple_angmom(basis.UC,'s_Na','s_Cs','S');
basis = rmfield(basis,'UC'); % don't need uncoupled basis anymore

%% projectors onto singlet and triplet subspaces
basis.SC.ops.Proj_singlet = 1*diag(basis.SC.qnums.S==0);
basis.SC.ops.Proj_triplet = 1*diag(basis.SC.qnums.S==1);

% spin dipole-dipole interactions give rise to the FB resonance
% include electron-electron spin and electron-nuclear spin
basis.SC.ops.Vdd0 = (c.mu0*c.uB^2/(4*pi*c.hartree*c.abohr^3))...
    *(c.gs_Na*c.gs_Cs*dd_operator(basis.SC.ops,'s_Na','s_Cs'));

% basis.SC.ops.F_z = basis.SC.ops.i_Na_z + basis.SC.ops.i_Cs_z + basis.SC.ops.N_z + basis.SC.ops.S_z;

%% form hund's case a & b bases
[basis.b,basis.SC,basis.change.SC_b] = couple_angmom(basis.SC,'N','S','J');

% transformation to hund's case a
basis.aUC.qnums = build_basis({'eta','s_Na','s_Cs','J','S','Lambda'},...
    {1,c.s_Na,c.s_Cs,0:1,0:1,0},[0 0 0 1 1 0],'a');
del = (basis.aUC.qnums.S>0 & basis.aUC.qnums.J==0) ...
    | (basis.aUC.qnums.S==0 & basis.aUC.qnums.J>0);
basis.aUC.qnums(del,:) = [];
basis.aUC.ops = build_operators(basis.aUC.qnums);
basis.change.aUC_b = operator_matrix(@case_b2a_element,{basis.aUC.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});

% transform operators into case a
f = setdiff(fields(basis.b.ops),fields(basis.aUC.ops));
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = basis.change.aUC_b * basis.b.ops.(f{i}) * basis.change.aUC_b';
end
basis.change.SC_aUC = basis.change.SC_b * (basis.change.aUC_b');

%% truncate basis to specific mtot (=4 for feshbach state)
basis = truncate_basis(basis,@(ops) ops.N_sq,0);

%% trap
waist = 1064e-9; % meter
Power = 10e-3; % watt
I0 = 2*Power/(pi*waist^2); % W/m^2
trap_depth = (c.apol_Cs*I0)/(2*c.eps0*c.c)/c.hartree;
atrap = 2*((c.abohr)/waist).^2;
Vtrap =@(r) trap_depth*(1-exp(-atrap*r.^2));

%% total hamiltonian
W =@(r) herm(NaCsXPES(r).*basis.SC.ops.Proj_singlet ...
    + NaCsaPES(r).*basis.SC.ops.Proj_triplet + basis.SC.ops.Vdd0./(r.^3) ...
    + Vtrap(r).*basis.SC.ops.I);

%% call the solver
[E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,c.mu_nacs/c.me);

%% display results
diff(E_out)*c.hartree/c.h * 1e-6
% errstr(E_out/c.wavenum2hartree*(29.9792458),err_est/c.wavenum2hartree*(29.9792458))

disp(E_out/c.wavenum2hartree*(29.9792458))

figure(2);
clf;
set(gcf,'units','inches','position',[5 5 3 2.4]);
for i = 1:size(psi,3)
    subplot(size(psi,3),1,i)
    plot(r,sum(abs(psi(:,:,i)).^2,1),'linewidth',2),set(gca,'xscale','log');
    set(gca,'fontsize',8)
    xlim([min(r) max(r)])
    ylabel('\psi(R)')
    if i<size(psi,3)
        set(gca,'xticklabel',[]);
    end
end
xlabel('R (a_0)')

%% change basis for output
out.E = E_out;
out.nodes = nodes_out;
out.r = r;
out.qnums = basis.(save_basis).qnums;
out.ops = basis.(save_basis).ops;
out.psi = zeros(size(basis.(save_basis).qnums,1),numel(r),size(basis.SC.qnums,3));
for i = 1:size(psi,3)
    out.psi(:,:,i) = (basis.change.(['SC_' save_basis])')*psi(:,:,i);
end
save(['data/feshbach_' datestr(now,'YYmmDD_HHMMSS') '.mat'],'out')
end