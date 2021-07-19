
function basis = nacs_basis(Ntot,Jmax)

if nargin<1
    Ntot = 0;
end
if nargin<2
    Jmax = Ntot+1;
end

Nmin = 0;
Nmax = Ntot+1;

basis.atom.qnums = build_basis({'s1','s2','l1','l2'},{1/2,1/2,0,0:1},[1 1 1 1]);
[basis.LC,basis.atom,basis.change.atom_LC] = couple_angmom(basis.atom,'l1','l2','L');

basis.LC.qnums.Lambda = basis.LC.qnums.m_L;
basis.LC.qnums = rmcol(basis.LC.qnums,{'L','m_L','l1'});
basis = rmfield(basis,'atom');

% molecular rotation
basis_mol.qnums = build_basis({'Lambda','N'},{[0 -1 1],Nmin:Nmax},[1,0],'b');
basis_mol.qnums(basis_mol.qnums.N>0 & basis_mol.qnums.Lambda==0,:) = [];

basis.LC = join_basis(basis.LC,basis_mol);
basis.LC.ops = build_operators(basis.LC.qnums);

[basis.SC,basis.LC,basis.change.SC_LC] = couple_angmom(basis.LC,'s1','s2','S');
basis = balance_operators(basis,'LC','SC');

[basis.b,basis.SC,basis.change.SC_b] = couple_angmom(basis.SC,'N','S','J');
basis.change.b_LC = basis.change.SC_b' * basis.change.SC_LC;
basis = balance_operators(basis,'SC','b');

basis = balance_operators(basis,'b','LC');

[basis.aUC,basis.change.aUC_b] = basis_b2a(basis.b);


basis.aUC.qnums.label(basis.aUC.qnums.S==0 & basis.aUC.qnums.Lambda==0 & basis.aUC.qnums.l2==1) = 'A';
basis.aUC.qnums.label(basis.aUC.qnums.S==0 & basis.aUC.qnums.Lambda==0 & basis.aUC.qnums.l2==0) = 'X';
basis.aUC.qnums.label(basis.aUC.qnums.S==0 & abs(basis.aUC.qnums.Lambda)==1 & basis.aUC.qnums.l2==1) = 'B';
basis.aUC.qnums.label(basis.aUC.qnums.S==1 & abs(basis.aUC.qnums.Lambda)==1 & basis.aUC.qnums.l2==1) = 'b';
basis.aUC.qnums.label(basis.aUC.qnums.S==1 & basis.aUC.qnums.Lambda==0 & basis.aUC.qnums.l2==1) = 'c';
basis.aUC.qnums.label(basis.aUC.qnums.S==1 & basis.aUC.qnums.Lambda==0 & basis.aUC.qnums.l2==0) = 'a';

term_letters = 'SPDF';
basis.aUC.qnums.term = strcat(basis.aUC.qnums.label,num2str(2*basis.aUC.qnums.S+1),...
    term_letters(abs(basis.aUC.qnums.Lambda)+1)',num2str(abs(basis.aUC.qnums.Omega)));

basis = truncate_basis(basis,@(ops) ops.J_sq,[0,Jmax*(Jmax+1)],0.1);

end