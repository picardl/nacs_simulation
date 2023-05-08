
const = constants();

Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;


cbB_file = '../data/cB_230209_131958.mat';
% fb_file = '../data/fb_860G_aFC_220623_121035.mat';%'../data/fb_852G_aFC_210706_113524.mat';
X_file = '../data/X_vib_210706_163735.mat';
a_file = '../data/a_210706_160129.mat';
% c_file = '../data/c_vib_210706_162602.mat';
emp_file = '../data/empirical_cbB_230204_205052.mat';

cbB_data = load(cbB_file);

E = cbB_data.out.E + Eoffs;
nodes = cbB_data.out.nodes;
r = cbB_data.out.r;
qnums = cbB_data.out.qnums;
ops = cbB_data.out.ops;
psi = cbB_data.out.psi;
W = cbB_data.out.W;

fb_data = feshbach(860e-4,'aFC',1,const);%load(fb_file);
r_fb = fb_data.r;
psi_fb = fb_data.psi(:,:,1);
qnums_fb = fb_data.qnums;

X_data = load(X_file);
r_X = X_data.out.r;
psi_X = X_data.out.psi;

a_data = load(a_file);
r_a = a_data.out.r;
psi_a = a_data.out.psi;

emp_data = load(emp_file);
nodes_c = emp_data.out_exp.nodes_c;
E_c = emp_data.out_exp.E_c;
r_c = emp_data.out_exp.r;
psi_c = emp_data.out_exp.psi_c;
nodes_B = emp_data.out_exp.nodes_B;
E_B = emp_data.out_exp.E_B;
r_B = emp_data.out_exp.r;
psi_B = emp_data.out_exp.psi_B;

E_emp = {E_c,E_B};

[V,W_adiab] = eigenshuffle(W);
Rref = 8;
ind_ref = find(min(abs(r-Rref))==abs(r-Rref));
[~,ind] = max(abs(V(:,:,ind_ref)).^2,[],1);
W_adiab = W_adiab(ind,:);

Nchn = size(qnums,1);
Nstates = size(psi,3);

%% lifetime calculation
lower_terms = {'X1S','a3S'};
Nlower = numel(lower_terms);

upper_terms = cellstr(qnums.term(:,1:3));
Nupper = numel(upper_terms);

Nr = numel(r);

w3d2 = zeros(Nupper,Nr);
for i = 1:Nupper
    for j = 1:Nlower
        fn = ['lib/rosario_potentials/DM_' lower_terms{j} '_' upper_terms{i} '_Full'];
        if exist(fn,'file')
            cbB_data = load(fn);
            lower_fun_name = ['NaCs' lower_terms{j}(1) 'PES'];
            if strcmp(upper_terms{i}(1),'b')
                upper_fun_name = ['NaCsbbPES'];
            else
                upper_fun_name = ['NaCs' upper_terms{i}(1) 'PES'];
            end
            Vlower = feval(lower_fun_name,r);
            Vupper = feval(upper_fun_name,r);
            omega_diff = const.hartree*(Vupper-Vlower)/const.hbar;
            dipole = const.e*const.abohr*interp1(cbB_data(:,1),cbB_data(:,2),r,'linear','extrap');
            w3d2(i,:) = w3d2(i,:) + omega_diff.^3 .* dipole.^2;
        end
    end
end
w3d2 = sum(trapz(r,abs(psi).^2.*w3d2,2),1);
gamma = w3d2/(3*pi*const.eps0*const.hbar*const.c^3);
gamma = gamma(:);

Rtest = 8.8;
ind_test = find(abs(r-Rtest)==min(abs(r-Rtest)));
[~,ind] = max(abs(psi(:,ind_test,:)).^2,[],1);

p = reshape(trapz(r,abs(psi).^2,2),Nchn,Nstates);
% [~,ind] = max(p,[],1);

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

%% franck-condon
for k = 1:Nchn
    dip_fb{k} = zeros(Nupper,Nlower,size(psi_chn{k},3));
    dip_X{k} = zeros(Nupper,Nlower,size(psi_chn{k},3),size(psi_X,3));
    dip_a{k} = zeros(Nupper,Nlower,size(psi_chn{k},3),size(psi_a,3));
    for i = 1:Nupper
        for j = 1:Nlower
            s = (str2num(lower_terms{j}(2))-1)/2;
            fn = ['lib/rosario_potentials/DM_' lower_terms{j} '_' upper_terms{i} '_Full'];
            if exist(fn,'file')
                cbB_data = load(fn);
                
                dipole = const.e*const.abohr*interp1(cbB_data(:,1),cbB_data(:,2),r,'linear','extrap');
                
                psi_lower = sum(interp1(r_fb(:),permute(psi_fb(qnums_fb.S==s,:,:),[2 1 3]),r(:)),2);
                integrand = psi_chn{k}(i,:,:).*dipole.*psi_lower';
                dip_fb{k}(i,j,:) = abs(trapz(r,integrand,2));
                
                if s==0
                    psi_lower = sum(interp1(r_X(:),permute(psi_X,[2 1 3]),r(:)),2);
                    integrand = psi_chn{k}(i,:,:).*dipole.*permute(psi_lower,[2 1 4 3]);
                    integrand(isnan(integrand)) = 0;
                    dip_X{k}(i,j,:,:) = abs(trapz(r,integrand,2));
                    
                    dip_a{k}(i,j,:,:) = 0;
                else
                    psi_lower = sum(interp1(r_a(:),permute(psi_a,[2 1 3]),r(:)),2);
                    integrand = psi_chn{k}(i,:,:).*dipole.*permute(psi_lower,[2 1 4 3]);
                    integrand(isnan(integrand)) = 0;
                    dip_a{k}(i,j,:,:) = abs(trapz(r,integrand,2));
                    
                    dip_X{k}(i,j,:,:) = 0;
                end
            end
        end
    end
    dip_fb{k} = sum(sum(dip_fb{k},2),1);
    dip_X{k} = sum(sum(dip_X{k},2),1);
    dip_a{k} = sum(sum(dip_a{k},2),1);
    dip_fb{k} = squeeze(dip_fb{k});
    dip_X{k} = squeeze(dip_X{k});
    dip_a{k} = squeeze(dip_a{k});
end

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
% for i = 1:3
%     psi_plot = yscale*(psi_chn{i}*peak2peak(E)/numel(nodes)/4 + reshape(E_chn{i},1,1,[])) - yoffs;
%     for j = 1:size(psi_chn{i},3)
%         pp = plot(r,psi_plot(:,:,j));
%         arrayfun(@(p,w) set(p,'color',get(w,'color')),pp,Wplot)
%     end
% end
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

figure(4);
clf;
hold on;
box on;
for i = 1:Nchn
    plot(yscale*E_chn{i},1e-6*gamma_chn{i}/(2*pi),'o')
end
hold off;
xlabel('E (THz)')
ylabel('\Gamma/(2\pi) (MHz)')
legend('c^3\Sigma','B^1\Pi','b^3\Pi','location','sw')

figure(5);
clf;
for k = 1:Nchn
    sp(1) = subplot(3,3,k);
    plot((0:numel(dip_fb{k})-1),dip_fb{k}/(const.e*const.abohr))
    xlabel(['v_{' qnums{k,'term'} '}'])
    xlim([0,numel(dip_fb{k})-1])
    ylabel(['\langle v_{' qnums{k,'term'} '} | D | FB \rangle (e a_0)'])
    
    sp(2) = subplot(3,3,k+3);
    imagesc(0:size(dip_X{k},1)-1,0:size(dip_X{k},2)-1,dip_X{k}'/(const.e*const.abohr))
    xlim([0,numel(dip_fb{k})-1])
    xlabel(['v_{' qnums{k,'term'} '}'])
    ylabel(['v_{X}'])
    colorbar
    
    sp(3) = subplot(3,3,k+6);
    imagesc(0:size(dip_a{k},1)-1,0:size(dip_a{k},2)-1,log10(dip_a{k}'/(const.e*const.abohr)))
    xlim([0,numel(dip_fb{k})-1])
    xlabel(['v_{' qnums{k,'term'} '}'])
    ylabel(['v_{a}'])
    colorbar
    
    linkaxes(sp,'x')
end

% energy offset vs v=0, compare ab initio with empirical
figure(6);
clf;
hold on;
box on;
for i = 1:Nchn
    subplot(Nchn,1,i)
    stem((E_chn{i})*yscale,0:numel(E_chn{i})-1)
    hold on
    stem((E_emp{i})*yscale,0:numel(E_emp{i})-1)
    xlim([285,335])
    ylabel('v')
    xlabel('E-E_0 (THz)')
    legend('perturbed','empirical')
end
hold off;


%% plot fig 7
lc = lines(3);

rd = ceil(p - max(p)+0.001);
sp = cumsum(rd,2);
c3s = find(rd(1,:) ==1);
b1p = find(rd(2,:) ==1);

%manual tweaks
swap_ind = [];%[34,35,37:39,41:61];%[30,32,33,35,36,38,39]
swap_ind_bc = [];%[12];%[34,35,37:39,41:61];%[30,32,33,35,36,38,39]
c3s(swap_ind) = [];
b1p(swap_ind_bc) = [];

c3sE = E(c3s);
b1pE = E(b1p);
c3sv = sp(1,c3s);
b1pv = sp(2,b1p);


% purity = 3*var(p);
% padd = p + repmat([-1;0;1],1,138);

plotSettings = struct();
plotSettings.LineStyle = '-';

%%
xrange = [0.997*min(c3sE),max(c3sE)]*yscale;
figure(7)
clf;
subaxis(3,1,1)
scatter(c3sE*yscale,p(1,c3s),'MarkerFaceColor',min(lc(1,:) + 0.4,1))
hold on
scatter(c3sE*yscale,p(2,c3s),'MarkerFaceColor',min(lc(2,:) + 0.4,1))
xlim(xrange)
x1 = xline(320.010,'-','Label','v=22');
x1.LabelVerticalAlignment = 'middle';
x1.LabelHorizontalAlignment = 'center';
x2 = xline(325.129,'-.','Label','v=26');
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'center';
legend('c^3\Sigma','B^1\Pi','location','w')
ylabel('State admixture')
ax = gca;
ax.XTick = [];
ax.FontSize = 10;
title('c^3\Sigma')
xline(334.369,'color','red');

subaxis(3,1,2)
scatter(b1pE*yscale,p(1,b1p),'MarkerFaceColor',min(lc(1,:) + 0.4,1))
hold on
scatter(b1pE*yscale,p(2,b1p),'MarkerFaceColor',min(lc(2,:) + 0.4,1))
xlim(xrange)
legend('c^3\Sigma','B^1\Pi','location','w')
ylabel('State admixture')
ax = gca;
ax.XTick = [];
title('B^1\Pi')
ax.FontSize = 10;
xline(334.369,'color','red');

linewidths = [39,27,70,12,10,12,17,42,120];
linewidths_err = [13,1,10,3,3,8,3,9,30];
linewidths_E = [288.699,306.497,309.260,310.625,311.986,313.341,320.002,323.866,325.121];
subaxis(3,1,3)
errorbar(linewidths_E,linewidths,linewidths_err,'Marker','o','MarkerFaceColor',min(lc(1,:) + 0.4,1),'LineStyle','none');
xlim(xrange);
ax = gca;
ylabel('Linewidth [MHz]')
x1 = xline(320.010,'-');
x2 = xline(325.129,'-.');
xlabel('Energy [THz]')
ax.FontSize = 10;