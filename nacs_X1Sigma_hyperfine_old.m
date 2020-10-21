function nacs_X1Sigma_hyperfine_old

%% parameters
i1 = 3/2; % na nuc spin
i2 = 7/2; % cs nuc spin
Nmax = 0; % max rotation state

Bv = 1.7396e9; % rot const
eQq1 = -0.097e6; % Na elec quadrupole const
eQq2 = 0.150e6; % Cs elec quadrupole const
c1 = 14.2;
c2 = 854.5;
c3 = 105.6;
c4 = 3941.8;
gr = 1e-3; % rotational g-factor
g1 = 1.478; % Na nuc g-factor
g2 = 0.738; % Cs nuc g-factor
sigma1 = 639.2e-6; % Na nuc spin shielding
sigma2 = 6278.7e-6; % Cs nuc spin shielding
muN = 762.2593285; % nuclear magneton
muE = 4.6*(3.354e-30); % NaCs ground-state molecule-frame dipole moment
alpha_p = 1872.1153; % parallel polarizability
alpha_s = 467.0375; % perpendicular polarizability
alpha0 = (2*alpha_s+alpha_p)/3;
delta_alpha = alpha_p - alpha_s;
T20_alpha = sqrt(2/3)*delta_alpha;
gamma_pol = 0; %acos(1/3)/2; % laser polarization ellipticity parameter

%% build operators in uncoupled basis
bases = {'UC','SC','FC','F1C','F2C'};

basis.UC.qnums = build_basis({'i1','i2','N'},{i1,i2,0:Nmax},[1 1 1]);
basis.UC.ops = build_operators(basis.UC.qnums);
Nstates = size(basis.UC.qnums,1);

%% successively couple angular momenta to form coupled bases
[basis.SC,basis.UC,basis_change.UC_SC] = couple_angmom(basis.UC,'i1','i2','I');
[basis.FC,basis.SC,basis_change.SC_FC] = couple_angmom(basis.SC,'I','N','F');
basis_change.UC_FC = basis_change.UC_SC*basis_change.SC_FC;
[basis.F1C,~,basis_change.UC_F1C] = couple_angmom(basis.UC,'i1','N','F1');
[basis.F2C,~,basis_change.UC_F2C] = couple_angmom(basis.UC,'i2','N','F2');

%% spin-spin scalar coupling, i1.i2
% do it in all bases to check consistency
for a = bases
    basis.(a{1}).ops.i1doti2 = op_dot(basis.(a{1}).ops,'i1','i2');
end

%% spin-rotation scalar coupling, i1.N and i2.N
basis.FC.ops.i1dotN = op_dot(basis.FC.ops,'i1','i2');
for a = bases
    basis.(a{1}).ops.i1dotN = op_dot(basis.(a{1}).ops,'i1','N');
    basis.(a{1}).ops.i2dotN = op_dot(basis.(a{1}).ops,'i2','N');
end

%% electric quadrupole coupling
j1_dot_j2_cpl =@(f,j1,j2) (f.*(f+1)-j1.*(j1+1)-j2.*(j2+1))/2;
EQ =@(F,N,I) -(3*j1_dot_j2_cpl(F,N,I).^2 + 3/2*j1_dot_j2_cpl(F,N,I) - I.*(I+1).*N.*(N+1))./(2*I.*(2*I-1).*(2*N-1).*(2*N+3));
basis.F1C.ops.EQ1 = diag(EQ(basis.F1C.qnums.F1,basis.F1C.qnums.N,basis.F1C.qnums.i1));
basis.F2C.ops.EQ2 = diag(EQ(basis.F2C.qnums.F2,basis.F2C.qnums.N,basis.F2C.qnums.i2));

%% electric dipole matrix elements
[row,col] = ndgrid(1:Nstates,1:Nstates);
row = row(:);
col = col(:);
wigD_element =@(Np,mNp,k,p,q,N,mN) (-1).^mNp.*sqrt((2*N+1).*(2*Np+1)).*w3j(Np,-mNp,k,p,N,mN).*w3j(Np,0,k,q,N,0);
dipole_mat=@(b,p,r,c) reshape(...
    wigD_element(b.qnums.N(r),b.qnums.mN(r),1,p,0,b.qnums.N(c),b.qnums.mN(c)).*krondelta_multi(b,{'N','mN'},r,c),...
    Nstates,Nstates);

% spherical components
basis.UC.ops.dipole_z = dipole_mat(basis.UC,0,row,col);
basis.UC.ops.dipole_p = dipole_mat(basis.UC,1,row,col);
basis.UC.ops.dipole_m = dipole_mat(basis.UC,-1,row,col);

% cartesian components
basis.UC.ops.dipole_x = (basis.UC.ops.dipole_m-basis.UC.ops.dipole_p)/sqrt(2);
basis.UC.ops.dipole_y = (basis.UC.ops.dipole_m+basis.UC.ops.dipole_p)*1i/sqrt(2);

%% convert all operators back to uncoupled basis for all subsequent calculations
UC_ops = fields(basis.UC.ops);
for a = 2:numel(bases)
    b_ops = fields(basis.(bases{a}).ops);
    not_in_UC = find(~ismember(b_ops,UC_ops))';
    for i = not_in_UC
        basis.UC.ops.(b_ops{i}) = basis_change.(['UC_' bases{a}])*basis.(bases{a}).ops.(b_ops{i})*basis_change.(['UC_' bases{a}])';
    end
end

%% light shift hamiltonian
% clear gamma_pol
% syms gamma_pol
% assume(gamma_pol,'real')
eps_odt = [cos(gamma_pol),1i*sin(gamma_pol),0];
T1_eps_odt = sphten(eps_odt);
T1_eps_odt_dag = sphten(conj(eps_odt));

T2_eps_odt = sphten_compose(2,T1_eps_odt_dag,T1_eps_odt);
T0_eps_odt = sphten_compose(0,T1_eps_odt_dag,T1_eps_odt);

ACstark_mat =@(b,r,c) reshape(krondelta_multi(b,{},r,c)*T0_eps_odt*alpha0 + ...
    sum(cell2mat(arrayfun(@(p) (-1).^p .* T2_eps_odt(3-p) .* wigD_element(b.qnums.N(r),b.qnums.mN(r),2,p,0,b.qnums.N(c),b.qnums.mN(c))*T20_alpha ,reshape(-2:2,1,1,[]),'un',0)),3).*krondelta_multi(b,{'N','mN'},r,c) ...
    ,Nstates,Nstates);
H_ACStark = ACstark_mat(basis.UC,row,col);
H_ACStark = (H_ACStark + H_ACStark')/2;

%% build hamiltonian in uncoupled basis
HUC =@(E,B,Elas) ...
    Bv * basis.UC.ops.Nsq + ...
    c4 * basis.UC.ops.i1doti2 + ...
    c1 * basis.UC.ops.i1dotN + ...
    c2 * basis.UC.ops.i2dotN + ...
    eQq1 * basis.UC.ops.EQ1 + ...
    eQq2 * basis.UC.ops.EQ2 + ...
    -gr*muN*( B * basis.UC.ops.Nz ) + ...
    -g1*muN*(1-sigma1)*( B * basis.UC.ops.i1z ) + ...
    -g2*muN*(1-sigma2)*( B * basis.UC.ops.i2z ) + ...
    -muE* ( E(1) * basis.UC.ops.dipole_x + E(2) * basis.UC.ops.dipole_y + E(3) * basis.UC.ops.dipole_z ) + ...
    + 0.5 * Elas^2 * H_ACStark;

%% tests
Isq_test = (basis.UC.ops.i1x+basis.UC.ops.i2x)^2 + (basis.UC.ops.i1y+basis.UC.ops.i2y)^2 + (basis.UC.ops.i1z+basis.UC.ops.i2z)^2;
test(1) = max(abs(reshape(basis_change.UC_SC'*Isq_test*basis_change.UC_SC - basis.SC.ops.Isq,[],1))) < 1e-10;
test(2) = max(abs(reshape(basis_change.UC_SC*basis.SC.ops.Isq*basis_change.UC_SC' - Isq_test,[],1))) < 1e-10;


EQiMatrixFC = basis_change.UC_FC'*(basis.UC.ops.EQ1 + basis.UC.ops.EQ2)*basis_change.UC_FC;

% [a,b] = evec_ind({'F','mF'},[3,3],basis.FC);

%% save a set of eigenvectors and a basis
B = 865;
Elas = 0;

H = cell2mat(reshape(arrayfun(@(x) HUC([0 0 0],x,Elas),B,'un',0),1,1,[]));

[V,D] = eigenshuffle(H);
D = real(D);

psi = V;
E = D;
qnums = basis.UC.qnums;
save(['NaCs_hfwfn_X1Sigma_' sprintf('%0.0f',B) 'G.mat'],'qnums','psi','B','E');

end

function j1dotj2 = op_dot(ops,j1_name,j2_name)

j1dotj2 = 0*ops.([j1_name 'x']);
for q = 'xyz'
    j1dotj2 = j1dotj2 + ops.([j1_name q])*ops.([j2_name q]);
end

end


