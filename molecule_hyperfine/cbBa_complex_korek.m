function out = cbBa_complex_korek()

const = constants();

Eoffs = -0.331957;

%% simulation parameters
Nx = 1e4;
rmin = 4.8;
rmax = 30;
% Erange = [-0.333 -0.331971]+0.01;
% Erange = [-0.2844 -0.2830] - Eoffs;
Erange = [324000 326000]*1e9*const.h/const.hartree;

%% basis
basis = nacs_basis(0,1);
basis = basis.aUC;

%% cut some terms
cut_terms = {'X1S0','c3S0','b3P2','a3S1','b3P0','A1S0','a3S0','a3S1'};
% cut_terms = {'X1S0','c3S0','b3P2','a3S1','a3S0','A1S0','b3P1'};
% cut_terms = {'X1S0','c3S0','b3P2','a3S1','A1S0','b3P1'};
basis.ops.cut = diag(ismember(strcat(basis.qnums.term,num2str(abs(basis.qnums.Omega))),cut_terms));
basis = truncate_basis(basis,@(ops) ops.cut,0,0.1);

basis.ops.Omega_z = diag(basis.qnums.Omega);
basis.ops.Lambda_z = diag(basis.qnums.Lambda);
basis = truncate_basis(basis,@(ops) ops.Omega_z,[0 1],0.1);
basis = truncate_basis(basis,@(ops) ops.J_z,0,0.1);
basis = truncate_basis(basis,@(ops) ops.J_sq - abs(ops.Omega_z).*(abs(ops.Omega_z)+ops.I),0,0.1);

%% load data
Nstates = size(basis.qnums,1);
R = linspace(rmin,rmax,Nx);

Wdiab = korek_potential(R,cellstr(basis.qnums.term));

Wint = zeros(Nstates,Nstates,numel(R));
for i = 1:Nstates
    Wint(i,i,:) = Wdiab(:,i);
end

terms = cellstr(basis.qnums.term(:,1:3));

NAC_A1S_B1P = load('../lib/rosario_potentials/NAC_A1S_B1P_Full');

[r,c] = ndgrid(1:Nstates,1:Nstates);
Omega = basis.qnums.Omega;
Lambda = basis.qnums.Lambda;

    function out = socpad(x)
        pp = csape(x(:,1),x(:,2),'variational');
        out = fnval(pp,R);
        out = reshape(out,1,1,numel(R));
    end

SOC_b3P_A1S = load('../lib/rosario_potentials/SOC_b3P_A1S_Full');
xi1 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'A1S'))...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'A1S')))...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==0)) ...
    .* socpad(SOC_b3P_A1S);

SOC_b3P_b3P = load('../lib/rosario_potentials/SOC_b3P_b3P_Full');
V0Pi = ((strcmp(terms(r),'b3P') & strcmp(terms(c),'b3P')) ...
    & (Omega(r)==-Omega(c)) & (Lambda(r)==-Lambda(c))...
    & (Omega(r)==0))  ...
    .* socpad(SOC_b3P_b3P);

SOC_b3P_B1P = load('../lib/rosario_potentials/SOC_b3P_B1P_Full');
zeta1 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'B1P'))) ...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==1)) ...
    .* socpad(SOC_b3P_B1P);

SOC_b3P_c3S = load('../lib/rosario_potentials/SOC_b3P_c3S_Full');
zeta2 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'c3S')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'c3S'))) ...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==1)) ...
    .* socpad(SOC_b3P_c3S);

SOC_c3S_B1P = load('../lib/rosario_potentials/SOC_c3S_B1P_Full');
zeta3 = (((strcmp(terms(r),'c3S') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'c3S') & strcmp(terms(r),'B1P'))) ...
    & ((Omega(r)==Omega(c)) & (abs(Omega(r))==1))) ...
    .* socpad(SOC_c3S_B1P);

dm_c3S_a3S = load('../lib/rosario_potentials/DM_a3S_c3S_Full');

Wint = Wint - zeta1 - zeta2 - zeta3 - V0Pi + xi1/sqrt(2) - Eoffs*eye(Nstates);

%%
[V,D] = eigenshuffle(Wint);
figure(1);
clf;
hold on;
plot(R*const.abohr*1e10,diag_nd(Wint),'--');
plot(R*const.abohr*1e10,D);
hold off;
set(gca,'xscale','log')

figure(3)
clf
semilogx(SOC_b3P_A1S(:,1)*const.abohr*1e10,SOC_b3P_A1S(:,2))
hold on
semilogx(SOC_b3P_b3P(:,1)*const.abohr*1e10,SOC_b3P_b3P(:,2))
semilogx(SOC_b3P_B1P(:,1)*const.abohr*1e10,SOC_b3P_B1P(:,2))
semilogx(SOC_b3P_c3S(:,1)*const.abohr*1e10,SOC_b3P_c3S(:,2))
semilogx(SOC_c3S_B1P(:,1)*const.abohr*1e10,SOC_c3S_B1P(:,2))
legend({'b3P-A1S','b3P-b3P','b3P-B1P','b3P-c3S','c3S-B1P'})

%% total hamiltonian
Wint = herm(Wint);
Wint = permute(reshape(Wint,Nstates^2,numel(R)),[2 1]);
W =@(r) permute(reshape(interp1(R,Wint,r,'linear','extrap'),1,1,numel(r),Nstates,Nstates),[4 5 3 1 2]);

%% call the solver
[E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,const.mu_nacs/const.me);

%% change basis for output
out.E = E_out;
out.nodes = nodes_out;
out.r = r;
out.W = W(reshape(r,1,1,[]));
out.qnums = basis.qnums;
out.ops = basis.ops;
out.psi = psi;
save(['../data/cbBa_' datestr(now,'YYmmDD_HHMMSS') '.mat'],'out')

end