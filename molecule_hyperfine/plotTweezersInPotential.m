
const = constants();

Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

cbB_file = '../data/cbB_210703_014725.mat';

X_file = '../data/X_vib_210706_163735.mat';
a_file = '../data/a_210706_160129.mat';
c_file = '../data/c_vib_210706_162602.mat';

cbB_data = load(cbB_file);

E = cbB_data.out.E + Eoffs;
nodes = cbB_data.out.nodes;
r = cbB_data.out.r;
qnums = cbB_data.out.qnums;
ops = cbB_data.out.ops;
psi = cbB_data.out.psi;
W = cbB_data.out.W;

X_data = load(X_file);
r_X = X_data.out.r;
psi_X = X_data.out.psi;

a_data = load(a_file);
r_a = a_data.out.r;
psi_a = a_data.out.psi;

c_data = load(c_file);
nodes_c = c_data.out.nodes;
E_c = c_data.out.E;
r_c = c_data.out.r;
psi_c = c_data.out.psi;

[V,W_adiab] = eigenshuffle(W);
Rref = 9;
ind_ref = find(min(abs(r-Rref))==abs(r-Rref));
[~,ind] = max(abs(V(:,:,ind_ref)).^2,[],1);
W_adiab = W_adiab(ind,:);

Nchn = size(qnums,1);
Nstates = size(psi,3);

E_chn = {};
psi_chn = {};
gamma_chn = {};
for i = 1:Nchn
    psi_chn{i} = psi(:,:,ind==i);
    E_chn{i} = E(ind==i);
    gamma_chn{i} = gamma(ind==i);
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

% figure(2);
% clf;
% hold on;
% box on;
% Wplot = plot(r,(W_adiab+Eoffs)*yscale-yoffs);
% for i = 1:3
%     psi_plot = yscale*(psi_chn{i}*peak2peak(E)/numel(nodes)/4 + reshape(E_chn{i},1,1,[])) - yoffs;
%     for j = 1:size(psi_chn{i},3)
%         pp = plot(r,psi_plot(:,:,j));
%         arrayfun(@(p,w) set(p,'color',get(w,'color')),pp,Wplot)
%     end
% end
% hold off;
% ylim([160 335]-yoffs)
% set(gca,'xscale','log')
% xlabel('R (a_0)')
% ylabel('E (THz)')
% set(gca,'fontsize',12)
% 
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
% legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','sw')
% 
% figure(5);
% clf;
% for k = 1:Nchn
%     sp(1) = subplot(3,3,k);
%     plot((0:numel(dip_fb{k})-1),dip_fb{k}/(const.e*const.abohr))
%     xlabel(['v_{' qnums{k,'term'} '}'])
%     xlim([0,numel(dip_fb{k})-1])
%     ylabel(['\langle v_{' qnums{k,'term'} '} | D | FB \rangle (e a_0)'])
%     
%     sp(2) = subplot(3,3,k+3);
%     imagesc(0:size(dip_X{k},1)-1,0:size(dip_X{k},2)-1,dip_X{k}'/(const.e*const.abohr))
%     xlim([0,numel(dip_fb{k})-1])
%     xlabel(['v_{' qnums{k,'term'} '}'])
%     ylabel(['v_{X}'])
%     colorbar
%     
%     sp(3) = subplot(3,3,k+6);
%     imagesc(0:size(dip_a{k},1)-1,0:size(dip_a{k},2)-1,log10(dip_a{k}'/(const.e*const.abohr)))
%     xlim([0,numel(dip_fb{k})-1])
%     xlabel(['v_{' qnums{k,'term'} '}'])
%     ylabel(['v_{a}'])
%     colorbar
%     
%     linkaxes(sp,'x')
% end
% 
% % energy offset vs v=0, compare ab initio with empirical
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
% 
% %% plot fig 7
% lc = lines(3);
% 
% rd = ceil(p - max(p)+0.001);
% sp = cumsum(rd,2);
% b3p = find(rd(3,:) ==1);
% c3s = find(rd(1,:) ==1);
% b1p = find(rd(2,:) ==1);
% 
% %manual tweaks
% % b3p_sup = c3s([30,32,33,35,36,38,39]);
% % c3s([30,32,33,35,36,38,39]) = [];
% % c3s_sup = b1p(12);
% % b1p(12) = [];
% 
% b3pE = E(b3p);
% c3sE = E(c3s);
% b1pE = E(b1p);
% b3pv = sp(3,b3p);
% c3sv = sp(1,c3s);
% b1pv = sp(2,b1p);
% 
% % b3pE_sup = E(b3p_sup);
% % c3sE_sup = E(c3s_sup);
% 
% 
% 
% % purity = 3*var(p);
% % padd = p + repmat([-1;0;1],1,138);
% 
% plotSettings = struct();
% plotSettings.LineStyle = '-';
% 
% xrange = [min(c3sE),max(c3sE)]*yscale;
% figure(7)
% clf;
% subaxis(3,1,3)
% scatter(b3pE*yscale,p(1,b3p),'MarkerFaceColor',min(lc(1,:) + 0.4,1))
% hold on
% scatter(b3pE*yscale,p(2,b3p),'MarkerFaceColor',min(lc(2,:) + 0.4,1))
% scatter(b3pE*yscale,p(3,b3p),'MarkerFaceColor',min(lc(3,:) + 0.4,1))
% 
% % scatter(b3pE_sup*yscale,p(1,b3p_sup),'MarkerEdgeColor',lc(1,:))
% % scatter(b3pE_sup*yscale,p(2,b3p_sup),'MarkerEdgeColor',lc(2,:))
% % scatter(b3pE_sup*yscale,p(3,b3p_sup),'MarkerEdgeColor',lc(3,:))
% xlim(xrange)
% legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','w')
% ylabel('State admixture')
% ax = gca;
% ax.FontSize = 10;
% xlabel('Energy [THz]')
% title('b^3\Pi')
% 
% subaxis(3,1,1)
% scatter(c3sE*yscale,p(1,c3s),'MarkerFaceColor',min(lc(1,:) + 0.4,1))
% hold on
% scatter(c3sE*yscale,p(2,c3s),'MarkerFaceColor',min(lc(2,:) + 0.4,1))
% scatter(c3sE*yscale,p(3,c3s),'MarkerFaceColor',min(lc(3,:) + 0.4,1))
% 
% % scatter(c3sE_sup*yscale,p(1,c3s_sup),'MarkerEdgeColor',lc(1,:))
% % scatter(c3sE_sup*yscale,p(2,c3s_sup),'MarkerEdgeColor',lc(2,:))
% % scatter(c3sE_sup*yscale,p(3,c3s_sup),'MarkerEdgeColor',lc(3,:))
% xlim(xrange)
% x1 = xline(318.757,'-','Label','v=22');
% x1.LabelVerticalAlignment = 'middle';
% x1.LabelHorizontalAlignment = 'center';
% x2 = xline(323.981,'-.','Label','v=26');
% x2.LabelVerticalAlignment = 'middle';
% x2.LabelHorizontalAlignment = 'center';
% legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','w')
% ylabel('State admixture')
% ax = gca;
% ax.XTick = [];
% ax.FontSize = 10;
% title('c^3\Sigma')
% 
% subaxis(3,1,2)
% scatter(b1pE*yscale,p(1,b1p),'MarkerFaceColor',min(lc(1,:) + 0.4,1))
% hold on
% scatter(b1pE*yscale,p(2,b1p),'MarkerFaceColor',min(lc(2,:) + 0.4,1))
% scatter(b1pE*yscale,p(3,b1p),'MarkerFaceColor',min(lc(3,:) + 0.4,1))
% xlim(xrange)
% legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','w')
% ylabel('State admixture')
% ax = gca;
% ax.XTick = [];
% title('B^1\Pi')
% ax.FontSize = 10;
