function plot_Ab_data

const = constants();

Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

Ab_file = 'data/Ab_210702_034846.mat';
fb_file = 'data/feshbach_210701_210331.mat';
X_file = 'data/X_210701_222231.mat';

Ab_data = load(Ab_file);

E = Ab_data.out.E + Eoffs;
nodes = Ab_data.out.nodes;
r = Ab_data.out.r;
qnums = Ab_data.out.qnums;
ops = Ab_data.out.ops;
psi = Ab_data.out.psi;
W = Ab_data.out.W;

fb_data = load(fb_file);
r_fb = fb_data.out.r;
psi_fb = fb_data.out.psi(:,:,1);
qnums_fb = fb_data.out.qnums;

X_data = load(X_file);
r_X = X_data.out.r;
psi_X = X_data.out.psi_r;
qnums_X = X_data.out.qnums;

[V,W_adiab] = eigenshuffle(W);
Rref = 7;
ind_ref = find(min(abs(r-Rref))==abs(r-Rref));
[~,ind] = max(abs(V(:,:,ind_ref)).^2,[],1);
W_adiab = W_adiab(ind,:);

psi = pagemtimes(V(:,:,ind_ref)',psi);

Nchn = size(qnums,1);
Nstates = size(psi,3);

%% lifetime calculation
lower_terms = {'X1S','a3S'};
Nlower = numel(lower_terms);

upper_terms = cellstr(qnums.term(:,1:3));
Nupper = numel(upper_terms);

Nr = numel(r);

% w3d2 = zeros(Nupper,Nr);
% for i = 1:Nupper
%     for j = 1:Nlower
%         fn = ['lib/rosario_potentials/DM_' lower_terms{j} '_' upper_terms{i} '_Full'];
%         if exist(fn,'file')
%             Ab_data = load(fn);
%             lower_fun_name = ['NaCs' lower_terms{j}(1) 'PES'];
%             if strcmp(upper_terms{i}(1),'b')
%                 upper_fun_name = ['NaCsbbPES'];
%             else
%                 upper_fun_name = ['NaCs' upper_terms{i}(1) 'PES'];
%             end
%             Vlower = feval(lower_fun_name,r);
%             Vupper = feval(upper_fun_name,r);
%             omega_diff = const.hartree*(Vupper-Vlower)/const.hbar;
%             dipole = const.e*const.abohr*interp1(Ab_data(:,1),Ab_data(:,2),r,'linear','extrap');
%             w3d2(i,:) = w3d2(i,:) + omega_diff.^3 .* dipole.^2;
%         end
%     end
% end
% w3d2 = sum(trapz(r,abs(psi).^2.*w3d2,2),1);
% gamma = w3d2/(3*pi*const.eps0*const.hbar*const.c^3);
% gamma = gamma(:);

p = reshape(trapz(r,abs(psi).^2,2),Nchn,Nstates);
[~,ind] = max(p,[],1);

E_chn = {};
psi_chn = {};
% gamma_chn = {};
for i = 1:Nchn
    psi_chn{i} = psi(:,:,ind==i);
    E_chn{i} = E(ind==i);
%     gamma_chn{i} = gamma(ind==i);
end

yscale = 1e-12*const.hartree/const.h;
yoffs = 0;

%%
figure(1);
clf;
hold on;
box on;
for i = 1:Nchn
    stem(E_chn{i}*yscale,0:numel(E_chn{i})-1)
end
hold off;
ylabel('v')
xlabel('E (THz)')
legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','sw')

figure(2);
clf;
hold on;
box on;
Wplot = plot(r,(W_adiab+Eoffs)*yscale-yoffs);
for i = 2
    psi_plot = yscale*(psi_chn{i}*peak2peak(E)/numel(nodes) + reshape(E_chn{i},1,1,[])) - yoffs;
    for j = 1:size(psi_chn{i},3)
        pp = plot(r,psi_plot(:,:,j));
        arrayfun(@(p,w) set(p,'color',get(w,'color')),pp,Wplot)
    end
end
hold off;
ylim([160 335]-yoffs)
set(gca,'xscale','log')
xlabel('R (a_0)')
ylabel('E (THz)')
set(gca,'fontsize',12)

figure(3);
clf;
stem(E*yscale,p','linewidth',1)
xlabel('E (THz)')
ylabel('State admixture')
legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','sw')

% figure(4);
% clf;
% hold on;
% box on;
% for i = 1:Nchn
%     plot(yscale*E_chn{i},1e-6*gamma_chn{i}/(2*pi),'o')
% end
% hold off;
% xlabel('E (THz)')
% ylabel('\Gamma/(2\pi) (MHz)')

end