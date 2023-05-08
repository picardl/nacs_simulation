function cB_deperturbation_lewis


const = constants();


colors = {'r','b','k'};

Eoffs = const.Cs_D12_weighted/const.hartree;

Nx = 1000;
xrange = [7.1 14];
Erange = [250 330]*1e12*const.h/const.hartree;
Eref = 320.010e12*const.h/const.hartree;

levelsssr = 0;

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
cut_terms = {'c3S0','b3P2','A1S0','b3P0','b3P1'};
% cut_terms = {'c3S0','b3P2'};
basis.qnums(ismember(basis.qnums.term,cut_terms),:) = [];
basis.qnums(basis.qnums.Sigma == -1,:) = [];
basis.ops = build_operators(basis.qnums);

basis.ops.Omega_z = diag(basis.qnums.Omega);
[basis,~,~] = truncate_basis(basis,@(ops) sign(ops.Omega_z),[1],1);

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


SOC_c3S_B1P = load('../lib/rosario_potentials/SOC_c3S_B1P_Full');
zeta3 = (((strcmp(terms(row),'c3S') & strcmp(terms(c),'B1P')) ...
    | (strcmp(terms(c),'c3S') & strcmp(terms(row),'B1P'))) ...
    & ((Omega(row)==Omega(c)) & (abs(Omega(row))==1))) ...
    .* socpad(SOC_c3S_B1P);

W_thy = Wint - zeta3;

% total hamiltonian
W_thy = reshape(W_thy,Nstates^2,numel(R));
    function out = W(r)
        pp = csape(R,W_thy,'variational');
        out = fnval(pp,r);
        out = reshape(out,Nstates,Nstates,numel(r));
    end

r = adaptive_grid(@W,const.m_nacs,Erange,xrange,Nx);
W_thy = W(r) + Eoffs*eye(Nstates);

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
        
        Mpert_temp = zeros(3,1,numel(r));
        for ii=1:3
            Mpert_temp(ii,1,:) = sum(Mpert(1,chn==ii,:),2);
        end
        Mpert = Mpert_temp;
        
        out = [
            Mpert(1,:,:), Mpert(3,:,:);
            Mpert(3,:,:),  Mpert(2,:,:);
            ];

    end

R_test = 7.5;
ind_test = abs(r-R_test)==min(abs(r-R_test));
W_thy(2,1,:) = W_thy(2,1,:) +0.000859;
W_thy(1,2,:) = W_thy(1,2,:) +0.000859;
W_thy(1,1,:) =W_thy(1,1,:) + 8.4362e-04;
W_thy(2,2,:) =W_thy(2,2,:) + 8.4362e-04;
[V_thy,W_adiab_thy] = eigenshuffle(W_thy);
[~,ind_max] = max(abs(V_thy(:,:,ind_test)).^2,[],1);
W_adiab_thy = W_adiab_thy(ind_max,:);

% R_test2 = 8.35;
% ind_test2 = abs(r-R_test2)==min(abs(r-R_test2));
% [V_thy2,W_adiab_thy2] = eigenshuffle(W_thy);
% [~,ind_max2] = max(abs(V_thy2(:,:,ind_test2)).^2,[],1);
% W_adiab_thy2 = W_adiab_thy2(ind_max2,:);


% deltaW = W_tweak(r);

W_exp = W_emp(r);
% [V,D] = eigenshuffle(W_exp);
% [~,ind_max_exp] = max(abs(V(:,:,ind_test)).^2,[],1);
% W_exp = W_exp(ind_max_exp,:,:);
W_adiab_exp = diag_nd(W_exp);
% figure(2);
% plot(r,weight(r));

% figure(6);
% clf;
% hold on;
% box on;
% for i = 1
%     stem((E_chn{i})*yscale,0:numel(E_chn{i})-1)
% end
% stem((E_c)*yscale,0:numel(E_c)-1)
% hold off;
% ylabel('v')
% xlabel('E-E_0 (THz)')
% legend('ab initio','empirical')

ssrLog = [];
 %% optimization
iter = 0;
    function ssr = ssrfun(X)        
        Wp = W_pert(r,X,X_cB);
        
        W_thy_i = W_thy;
        
        W_perturbed = W_thy_i + Wp;
        W_pert_diab = diag_nd(W_perturbed);
        W_perturbed(isinf(W_perturbed)) = 0;
        W_perturbed(isnan(W_perturbed)) = 0;
        
        
        [V_thy_new,W_thy_new] = eigenshuffle(W_perturbed);
        [~,ind_max] = max(abs(V_thy_new(:,:,ind_test)).^2,[],1);
        W_thy_new = W_thy_new(ind_max,:);
        resid = (W_adiab_exp - W_thy_new);
        ssr = 1e7*trapz(r,sum(resid.^2,1),2);
        
%         [~,ind_max2] = max(abs(V_thy_new(:,:,ind_test2)).^2,[],1);
%         W_thy_new2 = W_thy_new(ind_max2,:);
%         resid2 = (W_adiab_exp - W_thy_new2);
%         ssr2 = 1e7*trapz(r,sum(resid2.^2,1),2);
%         
%         ssr = ssr1 + ssr2;


        if ~mod(iter,10) || iter == 0
            figure(1);
            clf;
            sp(1) = subplot(2,1,1);
            hold on
            box on
            for ii = 1:2
                plot(r,W_thy_new(ii,:),['--' colors{ii}])
                plot(r,W_adiab_thy(ii,:),['-.' colors{ii}])
%                 plot(r,W_pert_diab(ii,:),['-.' colors{ii}])
                plot(r,W_adiab_exp(ii,:),['-' colors{ii}])
            end
            
            hold off
            ylim([-0.02 0.01]+Eoffs)
            set(gca,'xscale','log')
            xlabel('R (a_0)')
            ylabel('W(R) (E_h)')
            legend({'1','1 thy','1 emp','2','2 thy','2 emp','3','3 emp'})
            
            sp(2) = subplot(2,1,2);
            hold on;
            box on;
            for ii = 1:2
                plot(r,resid(ii,:),colors{ii});
            end
            hold off
            set(gca,'xscale','log')
            xlabel('R (a_0)')
            ylabel('residual (E_h)')
            legend({'1','2','3'})
            
            ylim([-1 1]*1e-3)
            
%             sp(3) = subplot(3,1,3);
%             hold on;
%             box on;
%             for ii = 1:3
%                 plot(r,resid2(ii,:),colors{ii});
%             end
%             hold off
%             set(gca,'xscale','log')
%             xlabel('R (a_0)')
%             ylabel('residual (E_h)')
%             legend({'1','2','3'})
%             
%             ylim([-1 1]*1e-3)
            
            figure(101);
            clf;
            for ii = 1:2
                for jj = 1:2
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
            
%             %             W2 = reshape(W_thy_new,Nstates^2,numel(r));
%             W = @(x) reshape(fnval(csape(r,W_perturbed,'variational'),x),Nstates,Nstates,numel(x));
%             [E_out_thy,~,~,~] = cc_logderiv_adaptive_multi(xrange,Nx,W,Erange,const.mu_nacs/const.me);
%             
%             %             nodesMin = min(length(E_out_thy),length(E_out_exp));
%             %             levelsssr = 2*sum(abs(E_out_exp(1:nodesMin) - E_out_thy(1:nodesMin)))/(peak2peak(Erange)/(max(max(W_thy_new)) - min(min(W_thy_new))));
%             if length(E_out_thy) > 0
%                 levelsssr = (E_out_thy(1) - Eref)/(peak2peak(Erange)/(max(max(W_thy_new)) - min(min(W_thy_new))));
%             else
%                 levelsssr = (peak2peak(Erange))/(peak2peak(Erange)/(max(max(W_thy_new)) - min(min(W_thy_new))));
%             end

            linkaxes(sp,'x')
            drawnow();          
            disp(ssr);
%             ssrLog(end+1) = ssr;
%             figure(3)
%             semilogy((1:length(ssrLog))*10,ssrLog);
%             ylabel('ssr')
%             xlabel('Iter')
        end
        
        ssr = ssr + levelsssr;
%         disp(ssr)
        iter = iter+1;
    end



% data = load('../data/cB_deperturbation_230205_200144.mat');
% data = load('../data/cB_deperturbation_230206_153632.mat');
data = load('../data/cB_deperturbation_230209_115051.mat');

X = reshape(data.X,4,[]);
X_cB = data.X_cB;
% X = X(:,1:4);
% X_cB = X_cB(1,1:4);
% % X = [X,X*0.1,X*0.01];
% % X_cB = [X_cB,X_cB,X_cB];
% X(:,7:8) = X(:,1:2) + X(:,1:2).*(randn(size(X(:,1:2)))/100);
% X_cB = [X_cB,1,2];
% % X(3:4,6) = X(3:4,6) + 5;
% % X(3:4,9) = X(3:4,9) + 5;

% X = X(:,[6,2,3,4]);
% X_cB = [1,2,3,3];
% X_cB = X_cB(1:6);
% X = X(:,[1,3,4,7,8]);
% X(:,5) = X(:,4);
% X(1,2:3) = -0.01*X(1,2:3);
% X_cB = X_cB([1,3,4,7,8]);

% X(1,3:6) = 0.1*X(1,3:6);
% X(2,[4,6]) = 7;
% X(3,3) = 4;
% X(1,1) = 20*X(1,1);
% X = X(:,[8,7,3,4]);
% X_cB = [1,2,3,3];
% X(1,1) = -0*X(1,1);
% X(1,2) = -0*X(1,2);
% X2 = zeros(4,4);
% X2(:,1) = X(:,7);
% X2(:,2) = [-0.0001,7,4,4];
% X(:,5:8) =  X(:,1:4);
% X(1,5:8) =  0.01*X(1,5:8);
% X = X + randn(size(X))/100;
% X = zeros(4,4);
% X(1,5:8) = randn(1,4)/1000-0.001;
% X(2,5:8) = randn(1,4) + 7;
% X(3:4,5:8) = randn(2,4)/2 + 2;
% X_cB = [1,2,3,3];
% X(1,2) = -0.001;
% X(3:4,2) =0.5;
% X(:,5) = [0.01,8.5,4,4];
% X(:,6) = [-0.005,7.2,5,5];
% X(1,:) = 0;
% X(1,1) = -0.01;
% X(2,1) = 9;
% X(4,1) = 3;
% X(3,1) = 2;
% X(1,2) = -0.01;
% X(2,2) = 7;
% X(4,2) = 2;
% X(3,2) = 2;
% X(:,5) = [-0.00001,5,1,1];
% X(:,6) = [-0.001,8,0.01,0.01];
% X(:,1) = [-0.001,8,1,1];
% X(:,2) = [-0.001,10,1,1];
% X_cB = [1,2,3,3];
% X(1,4) = -0.001;
% X_cB = [1,3,3,3];
% X = X2;

% X = X(:,X_cB==3);
% X_cB = X_cB(X_cB==3);
% X(1,2) = 0.05;
% X(4,2) = 5;
% X(2,3) = 6;
% X(1,3) = 0.01;
% X(3,3) = 3;
% X(4,3) = 3;
% % X(2,9) = 7;
% % X(1,[3,6,9]) = 0.001;
% % X = X(:,1:3);
% % X(:,4) = [-0.01;7.5;5;5];
% % % X(:,8) = 0.01*X(:,1);
% % X_cB = [1,2,3,3];
% % X(1,1) = 0.01;
% % X(2,1) = 7.5;
% % X(4,1) = 7;
% X(:,end+1) = X(:,2);
% X(1,end) = -0.001;
% X(2,end) = 8;
% X(3,end) = 10;
% X(4,end) = 10;
% X(:,end+1) = X(:,1);
% X(1,end) = -0.01;
% X(2,end) = 9;
% X(3,end) = 10;
% X(4,end) = 10;
% 
% X(:,end+1) = X(:,2);
% X(1,end) = 0.001;
% X(2,end) = 100;
% X(3,end) = 0.001;
% X(4,end) = 0.001;
% X(:,end+1) = X(:,1);
% X(1,end) = 0.001;
% X(2,end) = 100;
% X(3,end) = 0.01;
% X(4,end) = 0.01;
% 
% X_cB(end+1:end+4) = [2,1,2,1];





lb = ones(size(X))*-inf;
ub = ones(size(X))*inf;
lb(2:4,:) = 0;
% ub(3:4,1:2) = 0;
% ub(1,X_cB == 3) = abs(max(W_thy(1,2,:)));
% ub(1,1:2) = X(1,1:2);
% lb(1,1:2) = X(1,1:2);
% lb(4,:) = 0;

X = X(:)';

Nterms = size(X,2);

options = optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',5000);

% ssrfun(X(:)');
X = fminsearchbnd(@ssrfun,X(:)',lb(:)',ub(:)',options);
fn = ['../data/cB_deperturbation_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
save(fn,'X','X_cB','Wpert')
disp(fn)

Wpert_int = W_pert(r,X',X_cB);
Wint = W_thy + Wpert_int;

function fout = Wfun(r)
        pp = csape(r,Wint,'variational');
        fout = fnval(pp,r);
        fout = reshape(fout,Nstates,Nstates,numel(r));
end
%% call the solver
% [E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi(xrange,Nx,@Wfun,Erange,const.mu_nacs/const.me);
% out.E = E_out;
% out.nodes = nodes_out;
% out.r = r;
% out.W = W(reshape(r,1,1,[]));
% out.qnums = basis.qnums;
% out.ops = basis.ops;
% out.psi = psi;
% % 
% fn = ['../data/full_cB_deperturbation_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
% save(fn,'X','X_cB','Wpert','out')
% disp(fn)

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
        
        out = zeros(2,2,numel(r));
        
        out(1,1,:) = NaCscPES(r);
        out(2,2,:) = NaCsBPES(r);
    end


end