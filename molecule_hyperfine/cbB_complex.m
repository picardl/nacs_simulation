
function cbB_complex

const = constants();

% Coupled-channel calculation of excited cBb complex of NaCs

%% simulation parameters
Nx = 1001;
rmin = 5;
rmax = 1e2;
Erange = [-0.01 -0.003];

%% basis
basis.qnums = build_basis({'J','S','Lambda'},{0:2,0:1,0:1},[0 0 0],'a');
basis.qnums(basis.qnums.J ~= basis.qnums.m_J,:) = [];
basis.qnums(basis.qnums.J > abs(basis.qnums.Omega),:) = [];
basis.qnums = unique(basis.qnums,'rows');
[~,order] = sortrows(abs(basis.qnums{:,{'Omega','Lambda','S'}}));
basis.qnums = basis.qnums(order,:);

% term symbols
basis.qnums.label(basis.qnums.S==0 & basis.qnums.Lambda==0) = 'A';
basis.qnums.label(basis.qnums.S==0 & abs(basis.qnums.Lambda)==1) = 'B';
basis.qnums.label(basis.qnums.S==1 & abs(basis.qnums.Lambda)==1) = 'b';
basis.qnums.label(basis.qnums.S==1 & basis.qnums.Lambda==0) = 'c';

term_letters = 'SPDF';
basis.qnums.term = strcat(basis.qnums.label,num2str(2*basis.qnums.S+1),...
    term_letters(abs(basis.qnums.Lambda)+1)',num2str(abs(basis.qnums.Omega)));

%% cut some terms
% cut_terms = {'c3S0','b3P2','A1S0','b3P0','b3P1'};
cut_terms = {'c3S0','b3P2','A1S0','b3P0'};
% cut_terms = {'c3S0','b3P2'};
basis.qnums(ismember(basis.qnums.term,cut_terms),:) = [];
basis.ops = build_operators(basis.qnums);

basis.ops.Omega_z = diag(basis.qnums.Omega);
[basis,~,Nchn] = truncate_basis(basis,@(ops) sign(ops.Omega_z),1,1);

%% load data
terms = cellstr(basis.qnums.term(:,1:3));
uterms = unique(terms);
Nstates = size(basis.qnums,1);

term_ind =@(s) find(contains(terms,s));

data = load(['lib/rosario_potentials/A1S_FINAL_D']);
R = round(data(:,1),5);
R = R(:)';

Wint = zeros(Nstates,Nstates,numel(R));
for i = 1:numel(uterms)
    data = load(['lib/rosario_potentials/' uterms{i} '_FINAL_D']);
    term_ind_i = term_ind(uterms{i});
    for j = 1:numel(term_ind_i)
        Wint(term_ind_i(j),term_ind_i(j),:) = reshape(data(:,2),1,1,[]);
    end
end

NAC_A1S_B1P = load('lib/rosario_potentials/NAC_A1S_B1P_Full');
NRpad = sum(~ismember(R,round(NAC_A1S_B1P(:,1),5)));

[r,c] = ndgrid(1:Nstates,1:Nstates);
Omega = basis.qnums.Omega;
Lambda = basis.qnums.Lambda;

socpad = @(x) reshape(padarray(x(:,2),NRpad,0,'pre'),1,1,[]);


SOC_b3P_A1S = load('lib/rosario_potentials/SOC_b3P_A1S_Full');
xi1 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'A1S'))...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'A1S')))...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==0)) ...
    .* socpad(SOC_b3P_A1S);

SOC_b3P_b3P = load('lib/rosario_potentials/SOC_b3P_b3P_Full');
V0Pi = ((strcmp(terms(r),'b3P') & strcmp(terms(c),'b3P')) ...
    & (Omega(r)==-Omega(c)) & (Lambda(r)==-Lambda(c))...
     & (Omega(r)==0))  ...
    .* socpad(SOC_b3P_b3P);

SOC_b3P_B1P = load('lib/rosario_potentials/SOC_b3P_B1P_Full');
zeta1 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'B1P'))) ...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==1)) ...
    .* socpad(SOC_b3P_B1P);

SOC_b3P_c3S = load('lib/rosario_potentials/SOC_b3P_c3S_Full');
zeta2 = (((strcmp(terms(r),'b3P') & strcmp(terms(c),'c3S')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(r),'c3S'))) ...
    & (Omega(r)==Omega(c)) & (abs(Omega(r))==1)) ...
    .* socpad(SOC_b3P_c3S);

SOC_c3S_B1P = load('lib/rosario_potentials/SOC_c3S_B1P_Full');
zeta3 = (((strcmp(terms(r),'c3S') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'c3S') & strcmp(terms(r),'B1P'))) ...
    & ((Omega(r)==Omega(c)) & (abs(Omega(r))==1))) ...
    .* socpad(SOC_c3S_B1P);

Wint = Wint - zeta1 - zeta2 - zeta3 - V0Pi + xi1/sqrt(2);

[V,D] = eigenshuffle(Wint);
figure(1);
clf;
hold on;
% plot(R*const.abohr*1e10,diag_nd(Wint),'--');
plot(R*const.abohr*1e10,D);
hold off;
set(gca,'xscale','log')

%% total hamiltonian
Wint = herm(Wint);
Wint = reshape(Wint,Nstates^2,numel(R))';
W =@(r) permute(reshape(interp1(R,Wint,r),1,1,numel(r),Nstates,Nstates),[4 5 3 1 2]);

%% call the solver
% [E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,const.mu_nacs/const.me);

% %% change basis for output
% out.E = E_out;
% out.nodes = nodes_out;
% out.r = r;
% out.W = W(reshape(r,1,1,[]));
% out.qnums = basis.qnums;
% out.ops = basis.ops;
% out.psi = psi;
% fn = ['data/cbB_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
% save(fn,'out')
% disp(fn);

end