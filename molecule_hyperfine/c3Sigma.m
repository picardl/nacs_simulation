function out = c3Sigma(B,save_basis,recompute,Jmax,mtot,f_vib,c)

% Calculate rotational and hyperfine structure of c3Sigma state of NaCs


if nargin<1
    B = 860*1e-4; %B field in Tesla
end
if nargin<3
    recompute = 0;
end
if nargin<2
    save_basis = 'aFC';
end
if nargin<4
    Jmax = 3;
end
if nargin<5
    mtot = [3 4 5]; 
end
if nargin<6
    f_vib = 320010e9; %Target vibrational frequency in Hz
end
if nargin<7
    c = constants();
end

%% check for file at this B field with this basis
files = dir('../data');
file_ind = contains({files.name},['c_' strrep(num2str(B*1e4),'.','p') 'G_' save_basis]);
if any(file_ind) && ~(recompute > 1)
    disp('found c3Sigma file for this B field and basis')
    fnames = {files(file_ind).name};
    times = datenum(regexp(fnames,'\d{6}_\d{6}','match','once'),'YYmmDD_HHMMSS');
    data = load(['../data/' fnames{times==max(times)}]);
    out = data.out;
    return 
end

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

basis.aUC.ops.S_nz = operator_matrix(@spher_op,basis.aUC.qnums,{'S','Omega','Sigma'},0,'S','Omega');
basis.aUC.ops.S_np = operator_matrix(@spher_op,basis.aUC.qnums,{'S','Omega','Sigma'},1,'S','Omega');
basis.aUC.ops.S_nm = operator_matrix(@spher_op,basis.aUC.qnums,{'S','Omega','Sigma'},-1,'S','Omega');

basis.aUC.ops.Ss = 2*c.c3Sigma.lambda/3*(-basis.aUC.ops.S_np*basis.aUC.ops.S_np - basis.aUC.ops.S_nm*basis.aUC.ops.S_nm);

basis.aUC.ops.Hhf_Na = c.c3Sigma.alpha1*op_dot(basis.aUC.ops,'i_Na','S');
basis.aUC.ops.Hhf_Cs = c.c3Sigma.alpha2*op_dot(basis.aUC.ops,'i_Cs','S');

basis.aUC.ops.H = basis.aUC.ops.Hrot ...
    + basis.aUC.ops.HZ0elecspin.*B + basis.aUC.ops.Hhf_Na ...
    + basis.aUC.ops.Hhf_Cs + basis.aUC.ops.H_OmegaDoubling + basis.aUC.ops.Ss ;

%% go to fully coupled basis
[basis.aIC,basis.aUC,basis.change.aUC_IC] = couple_angmom(basis.aUC,'i_Na','i_Cs','I');
[basis.aFC,basis.aIC,basis.change.IC_FC] = couple_angmom(basis.aIC,'J','I','F');
basis.change.aUC_FC = basis.change.aUC_IC*basis.change.IC_FC;

%% transform F operators back to uncoupled basis
f = setdiff(fields(basis.aFC.ops),fields(basis.aUC.ops));
for i = 1:numel(f)
    basis.aUC.ops.(f{i}) = basis.change.aUC_FC*basis.aFC.ops.(f{i})*basis.change.aUC_FC';
end

%% truncate basis
% basis = rmfield(basis,'aUC');
% [basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.F_z,mtot);

%% diagonalize
[evecs,evals] = eig(basis.(save_basis).ops.H);
evals = real(diag(evals));

%% plot
mF_expect = real(diag(evecs'*basis.(save_basis).ops.F_z*evecs));
mJ_expect = real(diag(evecs'*basis.(save_basis).ops.J_z*evecs));
mI_expect = real(diag(evecs'*basis.(save_basis).ops.I_z*evecs));
figure(12);
clf;
hold on; box on;
for i = 1:numel(evals)
    plot(mF_expect(i) + [-0.4 0.4],[1 1]*evals(i)*1e-9/c.h,'-k')
end
hold off;
xlabel('M_{tot}');
ylabel('Energy (GHz)');

% figure(2);
% clf;
% hold on; box on;
% for i = 1:numel(evals)
%     plot(mI_expect(i) + [-0.4 0.4],[mJ_expect(i),mJ_expect(i)],'-k')
% end
% hold off;
% xlabel('M_{I}');
% ylabel('M_{J}');

%% load solution to deperturbed vib problem
vib_data = load('C:\nilab-projects\nacs_simulation\data\cbB_210703_014725.mat');
E_vib = vib_data.out.E*c.hartree + c.Cs_D12_weighted;

ind = find(abs(E_vib-f_vib*c.h) == min(abs(E_vib-f_vib*c.h)));

psi_vib = vib_data.out.psi(:,:,ind);
qnums_vib = vib_data.out.qnums;

%% save data
out.B = B;
out.qnums = basis.(save_basis).qnums;
out.ops = basis.(save_basis).ops;
out.psi = evecs;
out.r = vib_data.out.r;
out.nodes = vib_data.out.nodes(ind);
out.psi_vib = psi_vib;
out.qnums_vib = qnums_vib;
out.E = evals + E_vib(ind);

% fn = ['../data/c_' [strrep(num2str(B*1e4),'.','p') 'G'] '_vib' num2str(round(f_vib/1e12)) 'THz_' save_basis '_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
% save(fn,'out')
% disp(fn);

end
