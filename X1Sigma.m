
clear;


const = constants();

Nmax = 0;
gamma_pol = 0;

%% build operators in uncoupled basis
bases = {'UC','IC','FC','F1C','F2C'};

basis.UC.qnums = build_basis({'i_Na','i_Cs','N','S'},{const.i_Na,const.i_Cs,0:Nmax,0},[1 1 2 1]);
basis.UC.ops = build_operators(basis.UC.qnums);
Nstates = size(basis.UC.qnums,1);

%% successively couple angular momenta to form coupled bases
[basis.IC,basis.UC,basis_change.UC_IC] = couple_angmom(basis.UC,'i_Na','i_Cs','I');
[basis.FC,basis.IC,basis_change.IC_FC] = couple_angmom(basis.IC,'I','N','F');
basis_change.UC_FC = basis_change.UC_IC*basis_change.IC_FC;
[basis.F1C,~,basis_change.UC_F1C] = couple_angmom(basis.UC,'i_Na','N','F1');
[basis.F2C,~,basis_change.UC_F2C] = couple_angmom(basis.UC,'i_Cs','N','F2');
[basis.b,basis.UC,basis_change.UC_b] = couple_angmom(basis.UC,'N','S','J');

% transformation to hund's case a -- turns out it's the identity
basis.a.qnums = build_basis({'i_Na','i_Cs','J','S','Lambda'},{const.i_Na,const.i_Cs,0:Nmax,0,0},[1 1 2 0 0],'a');
basis.a.ops = struct();
basis_change.b_a = operator_matrix(@case_b2a_element,{basis.a.qnums,basis.b.qnums},{'J','Omega','S','Sigma','N','Lambda'});
basis_change.UC_a = basis_change.UC_b*basis_change.b_a;

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
        basis.UC.ops.(b_ops{i}) = basis_change.(['UC_' bases{a}])*basis.(bases{a}).ops.(b_ops{i})*basis_change.(['UC_' bases{a}])';
    end
end

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

B = 855e-4;

H = basis.UC.ops.H0 + basis.UC.ops.Hz0.*B;

[V,D] = eigenshuffle(H);
D = real(D);

