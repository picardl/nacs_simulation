function plot_cbB_data

const = constants();

Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

% cbB_file = 'data/cbB_210701_235922.mat';
% cbB_file = 'data/cbB_210703_014725.mat';
% cbB_file = 'data/cbB_210705_120344.mat';
% cbB_file = 'data/cbB_210705_125650.mat';
% cbB_file = 'data/cbB_210705_141528.mat';
% cbB_file = 'data/cbB_210706_154204.mat';
cbB_file = 'data/cbB_210708_103752.mat';

fb_file = 'data/fb_852G_aFC_210706_113524.mat';
X_file = 'data/X_vib_210706_163735.mat';
a_file = 'data/a_210706_160129.mat';
c_file = 'data/c_vib_210706_162602.mat';

cbB_data = load(cbB_file);

E = cbB_data.out.E + Eoffs;
nodes = cbB_data.out.nodes;
r = cbB_data.out.r;
qnums = cbB_data.out.qnums;
ops = cbB_data.out.ops;
psi = cbB_data.out.psi;
W = cbB_data.out.W;

fb_data = load(fb_file);
r_fb = fb_data.out.r;
psi_fb = fb_data.out.psi(:,:,1);
qnums_fb = fb_data.out.qnums;

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
for i = 1:3
    psi_plot = yscale*(psi_chn{i}*peak2peak(E)/numel(nodes)/4 + reshape(E_chn{i},1,1,[])) - yoffs;
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
for i = 1
    stem((E_chn{i})*yscale,0:numel(E_chn{i})-1)
end
stem((E_c)*yscale,0:numel(E_c)-1)
hold off;
ylabel('v')
xlabel('E-E_0 (THz)')
legend('ab initio','empirical')

end