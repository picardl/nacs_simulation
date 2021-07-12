function out = X1Sigma(B,save_basis,recompute,Nmax,mtot)

const = constants();

if nargin<1
    B = 853e-4;
end
if nargin<2
    recompute = 0;
end
if nargin<3
    save_basis = 'aFC';
end
if nargin<4
    Nmax = 2;
end
if nargin<5
    mtot = [2 3 4 5];
end

% vibrational solver parameters
Nx = 1e3;
rmin = 4; % abohr
rmax = 50; % abohr
Erange = -0.0224 + [-1 1]*1e-4; % energy range to search, atomic units
% Erange = [-0.0225 0];

%% check for file at this B field with this basis
files = dir('../data');
file_ind = contains({files.name},['X_' strrep(num2str(B*1e4),'.','p') 'G_' save_basis]);
if any(file_ind) && ~(recompute > 1)
    disp('found X1Sigma file for this B field and basis')
    fnames = {files(file_ind).name};
    times = datenum(regexp(fnames,'\d{6}_\d{6}','match','once'),'YYmmDD_HHMMSS');
    data = load(['../data/' fnames{times==max(times)}]);
    out = data.out;
    return 
end

%% build operators in uncoupled basis
basis.UC.qnums = build_basis({'eta','s_Na','s_Cs','Lambda','i_Na','i_Cs','N','S'},...
    {3,const.s_Na,const.s_Cs,0,const.i_Na,const.i_Cs,0:Nmax,0},[0 0 0 0 1 1 1 1]);
basis.UC.ops = build_operators(basis.UC.qnums);

%% successively couple angular momenta to form coupled bases
[basis.IC,basis.UC,basis.change.UC_IC] = couple_angmom(basis.UC,'i_Na','i_Cs','I');
[basis.FC,basis.IC,basis.change.IC_FC] = couple_angmom(basis.IC,'I','N','F');
[basis.F1C,basis.UC,basis.change.UC_F1C] = couple_angmom(basis.UC,'i_Na','N','F1');
[basis.F2C,basis.UC,basis.change.UC_F2C] = couple_angmom(basis.UC,'i_Cs','N','F2');
[basis.b,basis.UC,basis.change.UC_b] = couple_angmom(basis.UC,'N','S','J');

basis.change.UC_FC = basis.change.UC_IC*basis.change.IC_FC;
basis.change.F1C_FC = (basis.change.UC_F1C')*basis.change.UC_FC;
basis.change.F2C_FC = (basis.change.UC_F2C')*basis.change.UC_FC;

%% various operators
basis.UC.ops.F_z = basis.change.UC_FC * basis.FC.ops.F_z * basis.change.UC_FC';

% spin-spin scalar coupling, i_Na.i_Cs
basis.UC.ops.i_Nadoti_Cs = op_dot(basis.UC.ops,'i_Na','i_Cs');

% spin-rotation scalar coupling, i_Na.N and i_Cs.N
basis.UC.ops.i_NadotN = op_dot(basis.UC.ops,'i_Na','N');
basis.UC.ops.i_CsdotN = op_dot(basis.UC.ops,'i_Cs','N');

% electric quadrupole coupling
j1_dot_j2_cpl =@(f,j1,j2) (f.*(f+1)-j1.*(j1+1)-j2.*(j2+1))/2;
EQ =@(F,N,I) -(3*j1_dot_j2_cpl(F,N,I).^2 + 3/2*j1_dot_j2_cpl(F,N,I) - I.*(I+1).*N.*(N+1))./(2*I.*(2*I-1).*(2*N-1).*(2*N+3));
basis.F1C.ops.EQ1 = diag(EQ(basis.F1C.qnums.F1,basis.F1C.qnums.N,basis.F1C.qnums.i_Na));
basis.F2C.ops.EQ2 = diag(EQ(basis.F2C.qnums.F2,basis.F2C.qnums.N,basis.F2C.qnums.i_Cs));

% change to UC basis
basis.UC.ops.EQ1 = basis.change.UC_F1C * basis.F1C.ops.EQ1 * basis.change.UC_F1C';
basis.UC.ops.EQ2 = basis.change.UC_F1C * basis.F2C.ops.EQ2 * basis.change.UC_F1C';

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

basis.UC.ops.H = basis.UC.ops.H0 + basis.UC.ops.Hz0.*reshape(B,1,1,[]);

%% transformation to hund's case a
basis.aUC.qnums = build_basis({'eta','i_Na','i_Cs','J','S','Lambda'},{3,const.i_Na,const.i_Cs,0:Nmax,0,0},[0 1 1 2 0 0],'a');
basis.aUC.ops = struct();
basis.change.b_a = operator_matrix(@case_b2a_element,{basis.aUC.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});
basis.change.UC_a = basis.change.UC_b * basis.change.b_a;
f = fields(basis.UC.ops);
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = basis.change.UC_a'*basis.UC.ops.(f{i})*basis.change.UC_a;
end

[basis.aIC,basis.aUC,basis.change.a_aIC] = couple_angmom(basis.aUC,'i_Na','i_Cs','I');
[basis.aFC,basis.aIC,basis.change.aIC_aFC] = couple_angmom(basis.aIC,'J','I','F');

%% truncate basis
basis = rmfield(basis,{'IC','FC','F1C','F2C','b'});
[basis,~,~] = truncate_basis(basis,@(ops) ops.F_z,mtot);

%% solve for rotational/hyperfine energies
[psi_rot,E_rot] = eigenshuffle(basis.(save_basis).ops.H);
E_rot = real(E_rot);

% mF = diag(psi_rot'*basis.(save_basis).ops.F_z*psi_rot);
% figure(1);
% clf;
% hold on; box on;
% for i = 1:numel(E_rot)
%     plot([-1 1]*0.4+mF(i),[1 1]*E_rot(i)/const.h*1e-6,'-k');
% end
% hold off;
% xlabel('m_F');
% ylabel('E (MHz)');

%% solve for vibrational energy
[E_vib,nodes_out,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,@NaCsXPES,Erange,const.mu_nacs/const.me,1,1);

%% output
out.B = B;
out.r = r;
out.psi_r = psi_r;
out.nodes = nodes_out;
out.E = E_rot + E_vib*const.hartree;
out.psi = psi_rot;
out.qnums = basis.(save_basis).qnums;
out.ops = basis.(save_basis).ops;

fn = ['../data/X_' [strrep(num2str(B*1e4),'.','p') 'G'] '_' save_basis '_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
save(fn,'out')
disp(fn)

end
