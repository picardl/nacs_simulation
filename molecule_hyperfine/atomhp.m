function out = atomhp(B,save_basis,recompute,const)

if ~exist('c','var')
    const = constants();
end

if nargin<1
    B =8.8e-4;
end
if nargin<2
    save_basis = 'aFC';
end

% Coupled-channel calculation in ground state of NaCs.

if nargin<3
    recompute = 0;
end
if nargin < 4
    const = constants();
end
%% check for file at this B field with this basis
files = dir('../data');
file_ind = contains({files.name},['atomhp_' strrep(num2str(B*1e4),'.','p') 'G_' save_basis]);
if any(file_ind) && ~(recompute>1)
    disp('found atom hp channel file for this B field and basis')
    fnames = {files(file_ind).name};
    times = datenum(regexp(fnames,'\d{6}_\d{6}','match','once'),'YYmmDD_HHMMSS');
    data = load(['../data/' fnames{times==max(times)}]);
    out = data.out;
    return 
end

%% simulation parameters
Nx = 10000;
rmin = 4.5;
rmax = 10000;
N_tot = 0; % rotational quantum number
mtot = 6; % mtot = 6 for hp channel
Erange_vs_threshold = [-2e6 2e6]*const.h/const.hartree;%There should be no below threshold states of relevance

%% uncoupled basis
Nmin = 0;
Nmax = N_tot+1;

basis.UC = join_basis(atom_basis('Na'),atom_basis('Cs'));
basis_mol.qnums = build_basis({'eta','Lambda','N'},{1,0,Nmin:Nmax},[0,0,1],'b');
basis_mol.ops = build_operators(basis_mol.qnums);
basis.UC = join_basis(basis.UC,basis_mol);

% spin-coupled basis
[basis.SC,basis.UC,basis.change.UC_SC] = couple_angmom(basis.UC,'s_Na','s_Cs','S');
basis = rmfield(basis,'UC'); % don't need spin-uncoupled basis anymore

%% projectors onto singlet and triplet subspaces
basis.SC.ops.Proj_singlet = 1*diag(basis.SC.qnums.S==0);
basis.SC.ops.Proj_triplet = 1*diag(basis.SC.qnums.S==1);

% atomic hamiltonians
basis.SC.ops.H0 = (basis.SC.ops.H_Na_hyperfine + basis.SC.ops.H_Cs_hyperfine ...
    + (basis.SC.ops.H0_Na_zeeman + basis.SC.ops.H0_Cs_zeeman)*B)/const.hartree;

%% form hund's case b basis
basis.SC.ops.F_z = basis.SC.ops.i_Na_z + basis.SC.ops.i_Cs_z + basis.SC.ops.N_z + basis.SC.ops.S_z;
[basis.b,basis.SC,basis.change.SC_b] = couple_angmom(basis.SC,'N','S','J');

%% transformation to hund's case a
Jmax = Nmax; %max(basis.b.qnums.J);
Jmin = 0;
basis.aUC.qnums = build_basis({'eta','s_Na','s_Cs','i_Na','i_Cs','J','S','Lambda'},...
    {1,const.s_Na,const.s_Cs,const.i_Na,const.i_Cs,Jmin:Jmax,0:1,0},[0 0 0 1 1 0 0 0],'a');
basis.aUC.ops = build_operators(basis.aUC.qnums);
basis.change.aUC_b = operator_matrix(@case_b2a_element,{basis.aUC.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});

% transform operators into case a
f = setdiff(fields(basis.b.ops),fields(basis.aUC.ops));
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = pagemtimes(basis.change.aUC_b,pagemtimes( basis.b.ops.(f{i}),'none' , basis.change.aUC_b,'ctranspose'));
end
[basis.aIC,basis.aUC,basis.change.aUC_aIC] = couple_angmom(basis.aUC,'i_Na','i_Cs','I');
[basis.aFC,basis.aIC,basis.change.aIC_aFC] = couple_angmom(basis.aIC,'J','I','F');
basis.change.SC_aFC = basis.change.SC_b * (basis.change.aUC_b') * basis.change.aUC_aIC * basis.change.aIC_aFC;
basis.change.SC_aIC = basis.change.SC_b * (basis.change.aUC_b') * basis.change.aUC_aIC;
basis.change.SC_aUC = basis.change.SC_b * (basis.change.aUC_b');

%% truncate basis to specific mtot (=6 for hp channel)
basis = truncate_basis(basis,@(ops) ops.N_sq,N_tot*(N_tot+1));
basis = truncate_basis(basis,@(ops) ops.F_z,mtot);

%% rotation
basis.SC.ops.Hrot0 = const.hbar^2*basis.SC.ops.N_sq./(2*const.mu_nacs*const.hartree);

%% trap
waist = 1064e-9; % meter
Power = 10e-3;%10e-3; % watt
I0 = 2*Power/(pi*waist^2); % W/m^2
trap_depth = (const.apol_Cs*I0)/(2*const.eps0*const.c)/const.hartree;
atrap = 2*((const.abohr)/waist).^2;
Vtrap =@(r) trap_depth*(1-exp(-atrap*r.^2));

%% total hamiltonian
W =@(r) herm(basis.SC.ops.H0 + NaCsXPES(r).*basis.SC.ops.Proj_singlet ...
    + NaCsaPES(r).*basis.SC.ops.Proj_triplet + Vtrap(r).*basis.SC.ops.I + NaCsVd(r).*(basis.SC.ops.S_sq)) ;

%% find lowest energy channel and set energy search range
rtest = logspace(log10(rmin),log10(rmax),1e2);
Wtest = W(reshape(rtest,1,1,[]));
[V,D] = eigenshuffle(Wtest);
r_threshold = 300;
[E_lowest_chan_threshold,ind] = min(D(:,abs(rtest-r_threshold)==min(abs(rtest-r_threshold))));
Erange = E_lowest_chan_threshold + Erange_vs_threshold;

%% call the solver
[E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,const.mu_nacs/const.me);

%% change basis for output
out.B = B;
out.E = E_out*const.hartree;
out.nodes = nodes_out;
out.r = r;
out.qnums = basis.(save_basis).qnums;
out.ops = basis.(save_basis).ops;
out.psi = zeros(size(basis.(save_basis).qnums,1),numel(r),size(basis.SC.qnums,3));
out.E_lowest_chan_threshold = E_lowest_chan_threshold*const.hartree;
for i = 1:size(psi,3)
    out.psi(:,:,i) = (basis.change.(['SC_' save_basis])')*psi(:,:,i);
end

fn = ['../data/atomhp_' [strrep(num2str(B*1e4),'.','p') 'G'] '_' save_basis '_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
% save(fn,'out')
disp(fn)
end
