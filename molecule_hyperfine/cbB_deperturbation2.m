function cbB_deperturbation2


const = constants();


colors = {'r','b','k'};

Eoffs = const.Cs_D12_weighted/const.hartree;

Nx = 501;
xrange = [5 1e2];


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
[basis,~,~] = truncate_basis(basis,@(ops) sign(ops.Omega_z),1,1);

%% load data
terms = cellstr(basis.qnums.term(:,1:3));
uterms = unique(terms);
Nstates = size(basis.qnums,1);

term_ind =@(s) find(contains(terms,s));

data = load(['../lib/rosario_potentials/A1S_FINAL_D']);
R = round(data(:,1),5);
R = R(:)';

Wint = zeros(Nstates,Nstates,numel(R));
for i = 1:numel(uterms)
    data = load(['../lib/rosario_potentials/' uterms{i} '_FINAL_D']);
    term_ind_i = term_ind(uterms{i});
    for j = 1:numel(term_ind_i)
        Wint(term_ind_i(j),term_ind_i(j),:) = reshape(data(:,2),1,1,[]);
    end
end

NAC_A1S_B1P = load('../lib/rosario_potentials/NAC_A1S_B1P_Full');

[row,c] = ndgrid(1:Nstates,1:Nstates);
Omega = basis.qnums.Omega;
Lambda = basis.qnums.Lambda;

% socpad = @(x) reshape(padarray(x(:,2),NRpad,0,'pre'),1,1,[]);
    function out = socpad(x)
        pp = csape(x(:,1),x(:,2),'variational');
        out = fnval(pp,R);
        out = reshape(out,1,1,numel(R));
    end

SOC_b3P_A1S = load('../lib/rosario_potentials/SOC_b3P_A1S_Full');
xi1 = (((strcmp(terms(row),'b3P') & strcmp(terms(c),'A1S'))...
    | (strcmp(terms(c),'b3P') & strcmp(terms(row),'A1S')))...
    & (Omega(row)==Omega(c)) & (abs(Omega(row))==0)) ...
    .* socpad(SOC_b3P_A1S);

SOC_b3P_b3P = load('../lib/rosario_potentials/SOC_b3P_b3P_Full');
V0Pi = ((strcmp(terms(row),'b3P') & strcmp(terms(c),'b3P')) ...
    & (Omega(row)==-Omega(c)) & (Lambda(row)==-Lambda(c))...
     & (Omega(row)==0))  ...
    .* socpad(SOC_b3P_b3P);

SOC_b3P_B1P = load('../lib/rosario_potentials/SOC_b3P_B1P_Full');
zeta1 = (((strcmp(terms(row),'b3P') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(row),'B1P'))) ...
    & (Omega(row)==Omega(c)) & (abs(Omega(row))==1)) ...
    .* socpad(SOC_b3P_B1P);

SOC_b3P_c3S = load('../lib/rosario_potentials/SOC_b3P_c3S_Full');
zeta2 = (((strcmp(terms(row),'b3P') & strcmp(terms(c),'c3S')) ...
    | (strcmp(terms(c),'b3P') & strcmp(terms(row),'c3S'))) ...
    & (Omega(row)==Omega(c)) & (abs(Omega(row))==1)) ...
    .* socpad(SOC_b3P_c3S);

SOC_c3S_B1P = load('../lib/rosario_potentials/SOC_c3S_B1P_Full');
zeta3 = (((strcmp(terms(row),'c3S') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'c3S') & strcmp(terms(row),'B1P'))) ...
    & ((Omega(row)==Omega(c)) & (abs(Omega(row))==1))) ...
    .* socpad(SOC_c3S_B1P);

W_thy = Wint - zeta1 - zeta2 - zeta3 - V0Pi + xi1/sqrt(2);

% total hamiltonian
W_thy = reshape(W_thy,Nstates^2,numel(R));
    function out = W(r)
        pp = csape(R,W_thy,'variational');
        out = fnval(pp,r);
        out = reshape(out,Nstates,Nstates,numel(r));
    end

r = adaptive_grid(@W,const.m_nacs,xrange,Nx);
W_thy = W(r) + Eoffs*eye(3);

%%

Wpert =@W_pert;
    function out = W_pert(r,X,chn)
        X = reshape(X,4,[]);
        
        c = X(1,:);
        r0 = X(2,:);
        rpow = X(3,:);
        epow = X(4,:);
        
        r = reshape(r(:),1,1,1,[]);
        
        Mpert = real(c.*(r./(r0+eps)).^rpow.*exp(-(r./r0).^epow));
        Mpert = permute(Mpert,[3 2 4 1]);
        
        Mpert_temp = zeros(6,1,numel(r));
        for ii=1:6
            Mpert_temp(ii,1,:) = sum(Mpert(1,chn==ii,:),2);
        end
        Mpert = Mpert_temp;
        
        out = [
            Mpert(1,:,:) Mpert(4,:,:) Mpert(6,:,:);
            Mpert(4,:,:)  Mpert(2,:,:) Mpert(5,:,:);
            Mpert(6,:,:)   Mpert(5,:,:)   Mpert(3,:,:)
            ];
        
%         out(:,:,r>100) = 0;
    end

R_test = 8;
ind_test = abs(r-R_test)==min(abs(r-R_test));
[V_thy,W_adiab_thy] = eigenshuffle(W_thy);
[~,ind_max] = max(abs(V_thy(:,:,ind_test)).^2,[],1);
W_adiab_thy = W_adiab_thy(ind_max,:);

% deltaW = W_tweak(r);

W_exp = W_emp(r);
% [V,D] = eigenshuffle(W_exp);
% [~,ind_max_exp] = max(abs(V(:,:,ind_test)).^2,[],1);
% W_exp = W_exp(ind_max_exp,:,:);
W_adiab_exp = diag_nd(W_exp);

% figure(2);
% plot(r,weight(r));

%% optimization

iter = 0;
    function ssr = ssrfun(X)        
        Wp = W_pert(r,X,X_chn);
        
        W_thy_i = W_thy;
        
        W_perturbed = W_thy_i + Wp;
        W_pert_diab = diag_nd(W_perturbed);
        
        [V_thy_new,W_thy_new] = eigenshuffle(W_perturbed);
        [~,ind_max] = max(abs(V_thy_new(:,:,ind_test)).^2,[],1);
        W_thy_new = W_thy_new(ind_max,:);
        resid = (W_adiab_exp - W_thy_new);
        ssr = 1e7*trapz(r,sum(resid.^2,1),2);
        
        if ~mod(iter,100)
            figure(1);
            clf;
            sp(1) = subplot(2,1,1);
            hold on
            box on
            for ii = 1:3
                plot(r,W_thy_new(ii,:),['--' colors{ii}])
                plot(r,W_pert_diab(ii,:),['-.' colors{ii}])
                plot(r,W_adiab_exp(ii,:),['-' colors{ii}])
            end
            
            hold off
            ylim([-0.03 0.05]+Eoffs)
            set(gca,'xscale','log')
            xlabel('R (a_0)')
            ylabel('W(R) (E_h)')
            
            sp(2) = subplot(2,1,2);
            hold on;
            box on;
            for ii = 1:3
                plot(r,resid(ii,:),colors{ii});
            end
            hold off
            set(gca,'xscale','log')
            xlabel('R (a_0)')
            ylabel('residual (E_h)')
            legend('1','2','3')
            
            ylim([-1 1]*1e-3)
            
            figure(101);
            clf;
            for ii = 1:3
                for jj = 1:3
                    sp(1) = subplot(3,3,(ii-1)*3+jj);
                    hold on
                    box on
                    plot(r,squeeze(W_thy_i(ii,jj,:)),'-')
                    plot(r,squeeze(W_perturbed(ii,jj,:)),'-')
                    if ii==jj
                        plot(r,squeeze(W_exp(ii,jj,:)),'-')
                        ylim([-0.03 0.05]+Eoffs)
                    end
                    hold off
                                    
                    set(gca,'xscale','log')
                    xlabel('R (a_0)')
                    ylabel(['W_{' num2str(ii) num2str(jj) '}(R) (E_h)'])
                end
            end
            
            linkaxes(sp,'x')
            drawnow();
            
            disp(ssr)
            
            ssr;
        end
        
        
        iter = iter+1;
    end


% data = load('../data/deperturbation_210705_142800.mat');
% data = load('../data/deperturbation_210706_145339.mat');
% data = load('../data/deperturbation_210705_140759.mat');
% data = load('../data/deperturbation_210706_150901.mat');
% data = load('../data/deperturbation_210705_150222.mat');
data = load('../data/deperturbation_210712_141340.mat');

X = reshape(data.X,4,[]);
X_chn = data.X_chn;

% X = X(:,1:end-4);
% X_chn = X_chn(:,1:end-4);

% del = abs(X(4,:))>30;
% X(:,del) = [];
% X_chn(del) = [];

% c r0 rpow epow

% X = cat(2,X,[
%     1e-3 1e-3
%     7.2  7.2
%     10  10
%     10  10
%     ]);
% 
% X_chn = cat(2,X_chn,[4 1]);

lb = ones(size(X))*-inf;
ub = ones(size(X))*inf;
lb(2,:) = 0;
% ub(2,:) = max(r);
lb(4,:) = 0;

X = X(:)';

Nterms = size(X,2);

% ssrfun(X(:)');
X = fminsearchbnd(@ssrfun,X(:)',lb(:)',ub(:)');
% 
fn = ['../data/deperturbation_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
save(fn,'X','X_chn','Wpert')
disp(fn)

%%

% figure(1);
% clf;
% hold on
% box on
% for i = 1:3
%     plot(r,W_adiab_thy(i,:),['-' colors{i}])
%     plot(r,W_adiab_exp(i,:),['--' colors{i}])
% end
% hold off
% ylim([-0.03 0.005]+Eoffs)
% % set(gca,'xscale','log')

    function out = W_emp(r)
        r = reshape(r(:),1,1,[]);
        
        out = zeros(3,3,numel(r));
        
        out(1,1,:) = NaCscPES(r);
        out(2,2,:) = NaCsBPES(r);
        out(3,3,:) = NaCsbbPES(r);
    end


end