% [freqs,survivals,survivals_errs] =  loadFiles({'names_20210626_004412','names_20210626_163345','names_20210626_180552','names_20210626_234330','names_20210627_162921'},[1:66,115:129,137:215]);
% offs = 325125;
% freqs = freqs - offs;
% [freqs,order] = sort(freqs);
% survivals = survivals(order,:);
% survivals_errs = survivals_errs(order,:);

BFields = [10,433,860];
BField_labels = {'B = 10 G','B = 433 G','B = 860 G'};
% survivals = survivals_all{4}(1,:);
% survivals_errs = survivals_errs_all{4}(1,:);
offs = 325125;
% freqs = freqs_all{4} - offs;
r = constants;

r.B = 860*1e-4; %B Field in Gauss. updated FR in model is at 864.1
r.initState = 2; %1 for FB mol, 2 for atoms
r.basis = 'aFC';
r.recompute = 2;
r.Jmax = 3;%Max J to inlcude in rotational basis
r.mtot = [3:6]; %Range of mF to include in hyperfine basis
r.f_vib = 325118e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
r.pol_up = 'sigpm';
r.pol_frac = 0.5;
r.c3Sigma.Be = 1.04*1e9*r.h; %10.1103/PhysRevA.105.063322
r.c3Sigma.alpha1 =  0.7*1e9*r.h; %0.340*1e9*c.h;
r.c3Sigma.alpha2 =  0.09*1e9*r.h;
r.c3Sigma.Gamma = 100e6*r.h;
r.c3Sigma.gS = 2.0023;
r.power = 100e-3;
r.waist = 30e-6;
r.N_tot = 0;
vib_data = load('../data/cbB_210703_014725.mat');
E_vib = vib_data.out.E*r.hartree + r.Cs_D12_weighted;
ind = find(abs(E_vib-r.f_vib*r.h) == min(abs(E_vib-r.f_vib*r.h)));
r.E_vib = E_vib(ind);

f_off = 3;
load('../data/precomp_all_v.mat');
pre860 = prev26;
pre860.psi_init = zeros(size(pre860.psi_init));
pre860.psi_init(4) = 1;
fitfun860 = @(b,x) HamnRSC(x,r,pre860,pre860.f.E,b(1),b(2),b(3),b(4),b(5),r.pol_frac,mean(survivals(end-10:end,3))-min(survivals(:,3)),min(survivals(:,3)));

r.B = 433e-4;
pre433 = precomp_ops(r);
% pre433.psi_init = zeros(size(pre433.psi_init));
% pre433.psi_init(3) = 1;
% fitfun433 = @(b,x) HamnRSC(x,r,pre433,pre860.f.E,b(1),b(2),b(3),b(4),b(5),r.pol_frac,mean(survivals(end-10:end,2))-min(survivals(:,2)),min(survivals(:,2)));

r.B = 10e-4;
% pre10 = precomp_ops(r);
pre10.psi_init = zeros(size(pre10.psi_init));
pre10.psi_init(3) = 1;
fitfun10 = @(b,x) HamnRSC(x,r,pre10,pre860.f.E,b(1),b(2),b(3),b(4),b(5),r.pol_frac,mean(survivals(end-10:end,1))-min(survivals(:,1)),min(survivals(:,1)));

fitfuns = {fitfun10,fitfun433,fitfun860};
% 
%%
cols = lines(3);
figure(4)
clf
ftfreqs = linspace(min(freqs)-1,max(freqs),200);
hold on
i = 1;
    fitfun = fitfuns{i};
    ax = subaxis(3,2,i*2-1,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.08,'MarginRight',0.04,'MarginBottom',0.06,'MarginTop',0.03 );
    errorbar(freqs+125,survivals(:,i),survivals_errs(:,i),'Marker','o','markersize',5,'MarkerFaceColor',cols(1,:) + [0.1,0.1,0.1],'LineStyle','-','MarkerEdgeColor',cols(1,:),'Color',cols(1,:))
    xlim([124,140])
    legend(BField_labels(i*2-1),'location','nw')
    ax.XTickLabel = {};
    ylabel('Population')
    ylim([0.2,1])
    title('Data')
    ax.FontSize = 12;
    ax = subaxis(3,2,i*2,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.06,'MarginRight',0.02,'MarginBottom',0.06,'MarginTop',0.03  );
    plot(ftfreqs+125,fitfun([0.9511    0.3685   0.0212   -1.9675  0.0715],ftfreqs),'Color',cols(1,:),'LineWidth',2)
    xlim([124,140])
%     legend(BField_labels(i*2))
    ax.XTickLabel = {};
    ax.YTickLabel = {};
    ylim([0.2,1])
    title('Model')
    ax.FontSize = 12;
i = 2;
    fitfun = fitfuns{i};
ax = subaxis(3,2,i*2-1,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.08,'MarginRight',0.04,'MarginBottom',0.06,'MarginTop',0.01 );
errorbar(freqs+125,survivals(:,i),survivals_errs(:,i),'Marker','o','markersize',5,'MarkerFaceColor',cols(2,:)+ [0.1,0.2,0.2],'LineStyle','-','MarkerEdgeColor',cols(2,:),'Color',cols(2,:))
xlim([124,140])
legend(BField_labels(i*2-1),'location','nw')
ax.XTickLabel = {};
ylabel('Population')
ylim([0.2,1])
ax.FontSize = 12;
ax = subaxis(3,2,i*2,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.06,'MarginRight',0.02,'MarginBottom',0.06,'MarginTop',0.01  );
plot(ftfreqs+125,fitfun([0.9511    0.3685   0.0212   -1.9675  0.0715],ftfreqs),'Color',cols(2,:),'LineWidth',2)
xlim([124,140])
% legend(BField_labels(i*2))
ax.XTickLabel = {};
ax.YTickLabel = {};
ylim([0.2,1])
ax.FontSize = 12;
i = 3;
    fitfun = fitfuns{i};
    ax = subaxis(3,2,i*2-1,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.08,'MarginRight',0.04,'MarginBottom',0.06 );
    errorbar(freqs+125,survivals(:,i),survivals_errs(:,i),'Marker','o','markersize',5,'MarkerFaceColor',cols(3,:)+ [0.05,0.2,0.2],'LineStyle','-','MarkerEdgeColor',cols(3,:),'Color',cols(3,:))
    xlim([124,140])
    legend(BField_labels(i*2-1),'location','nw')
    ylabel('Population')
xlabel('325xxx [GHz]')
    ylim([0.2,1])
    ax.FontSize = 12;
    ax = subaxis(3,2,i*2,'Spacing', 0.02, 'Padding',  0.005, 'Margin', 0,'PaddingBottom',0.01,'PaddingTop',0,'MarginLeft',0.06,'MarginRight',0.02,'MarginBottom',0.06 );
    plot(ftfreqs+125,fitfun([0.9511    0.3685   0.0212   -1.9675  0.0715],ftfreqs),'Color',cols(3,:),'LineWidth',2)
    xlim([124,140])
%     legend(BField_labels(i*2))
    ax.YTickLabel = {};
    ylim([0.2,1])
xlabel('325xxx [GHz]')
ax.FontSize = 12;

%%
% figure(4)
% clf
% ftfreqs = linspace(min(freqs)-1,max(freqs),400);
% ftSurvs = zeros(400,3);
% hold on
% for i = 1:3
%     fitfun = fitfuns{i};
%     ftSurvs(:,i) = fitfun([0.9511    0.3685   0.0212   -1.9675  0.0715],ftfreqs);
% end
% subplot(1,2,1)
% stackedplot(freqs,survivals)
% xlabel('Freq [GHz]')
% xlim([-1,14])
% subplot(1,2,2)
% stackedplot(ftfreqs,ftSurvs)
% xlabel('Freq [GHz]')
% xlim([-1,14])

if 0
    %%
    r.c3Sigma.Be = beta(1)*1e9*r.h;
    r.c3Sigma.alpha1 = beta(2)*1e9*r.h;
    r.c3Sigma.alpha2 =  beta(3)*1e9*r.h;
    r.pol_up = 0.5;
    out = upleg_hamn_precomp(pre,r,r.B); 
    mJFreq = offToMj1(r,out,0.0910,offs)
end

if 0
    %%
    options = optimset;
    opstat = statset;
    opstat.Display = 'iter';
    options.Display = 'iter';
%     fitfun2 = @(b,x) fitfun([b(1),0.74,-0.11,b(2),b(3)],x);
%     lsqfitfun2 = @(b) lsqfitfun_t([b(1),b(2),b(3),b(4),b(5)]);

%     fitvals = fminsearch(lsqfitfun,[1   0.7105   -0.0253   0.6    0.5],options)
    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(freqs,survivals,fitfun,[0.95    0.6353   -0.5641    0.0830    0.004],opstat,'weights',1./survivals_errs.^2);
%     [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(freqs,survivals,fitfun2,[1.01   -1.1    0.0338],opstat,'weights',1./survivals_errs.^2);

end


function  popn = HamnRSC(f,r,pre,E_ref,Be,alpha1,alpha2,f_offs,t,pol_up,scale,y_offs)
    % %Compute effective Hamiltonian for up-leg (and down-leg) including
    % %rovibrational contribution to TDM

    r.c3Sigma.Be = Be*1e9*r.h;
    r.c3Sigma.alpha1 = alpha1*1e9*r.h;
    r.c3Sigma.alpha2 = alpha2*1e9*r.h;
    r.pol_up = pol_up;
    out = upleg_hamn_precomp(pre,r,r.B);
    popn = findRsc(f,out.E,pre.f.E,E_ref,out.rabi,r,t,f_offs,scale,y_offs,r.c3Sigma.Gamma);
end

function popn = findRsc(f,E,E_fb,E_ref,rabi,const,t,f_offs,scale,y_offs,Gamma)
    %% 
    %Calculate up-leg rabi rates to each state and corresponding total
    %scattering rates
    freq_offs = min(min(E - max(E_ref)))/const.h + f_offs*1e9;
    delta = E - E_fb - const.h*reshape(freq_offs + f*1e9,1,1,[]);
    Rsc = permute(sum(sum((Gamma/(2*const.h)) * 2*abs(rabi).^2./(2*abs(rabi).^2 + 4*delta.^2 + Gamma.^2),2),1),[2 3 1]);
    popn = abs(scale)*exp(-abs(Rsc*t)) + y_offs;
%     popn = exp(-Rsc*t);
end

function [freqs,survivals,survivals_errs] =  loadFiles(namesfile,inds)
    datadir = 'C:\NaCs1pt5_Data_temp\Data\';    
    
    
    if ~iscell(namesfile)
        try
            nameFiles = load([datadir namesfile(7:14) '\' namesfile '.mat']).names;
        catch
            nameFiles =  load([datadir, namesfile(7:14),'\', 'data_', namesfile,'.mat']); % load(FindDataFile(nameFiles{i})); % load([datadir{j},'data_', nameFiles{i},'.mat']);
        end
    else
        nameFiles = {};
        for i = 1:length(namesfile)
            try
                nameFiles = [nameFiles, load([datadir namesfile{i}(7:14) '\' namesfile{i} '.mat']).names];
            catch
                nameFiles = [nameFiles, load([datadir namesfile{i}(7:14) '\' namesfile{i} '.mat']).names];
            end
        end
    end 
    if ~exist('inds','var')
        inds = 1:length(nameFiles);
    end
    nameFiles = nameFiles(inds);
    
    N = size(nameFiles,2);
    freqs = zeros(1,N);
    survivals = zeros(1,N);
    survivals_errs = zeros(1,N)';
    
    for i = 1:N

        try
            loaddat =  load([datadir, nameFiles{i}(1:8),'\', 'data_', nameFiles{i},'\', 'data_', nameFiles{i},'.mat']);
        catch
            loaddat =  load([datadir, nameFiles{i}(1:8),'\', 'data_', nameFiles{i},'.mat']);
        end

        try
            survivals(1:7,i,1:3) = loaddat.Analysis.SurvivalProbability;
            survivals_errs(1:7,i,1:3) = loaddat.Analysis.SurvivalProbabilityErr;
        catch
            disp(['error loading file ' nameFiles{i}]);
            survivals(:,i,:) = nan;
            survivals_errs(:,i,:) = nan;
        end
        scgrp = ScanGroup.load(loaddat.Scan.ScanGroup);
        try
            freqs(i) = scgrp.get_fixed(1).nacs.pa.freq;
        catch
            if isfield(scgrp.get_fixed(1),'fWavemeter1')
                freqs(i) = scgrp.get_fixed(1).fWavemeter1;
            elseif isfield(scgrp.get_fixed(1),'fWavemeter')
                freqs(i) = (scgrp.get_fixed(1).fWavemeter) + 306000;
            end
        end
    end
    
    
    
    [freqs,order] = sort(freqs);
    survivals = squeeze(survivals(3,order,:));
    survivals_errs = squeeze(survivals_errs(3,order,:));
end

function pOut = polSwitch(pol)
    switch pol
        case 'sigp'
            pOut = [1 -1i 0]/sqrt(2);
        case 'sigm'
            pOut = [1 1i 0]/sqrt(2);
        case 'pi'
            pOut = [0 0 1];
        case 'sigpm'
            pOut = [1 0 0];
    end
end

function pOut = polFrac(sigp_frac)
    pOut = sqrt(sigp_frac)*[1 -1i 0]/sqrt(2) + sqrt(1 - sigp_frac)*[1 1i 0]/sqrt(2);
end

function mJFreq = offToMj1(r,out,fitOffs,fixedOffs)
    %fit offs is the output of the fit, fixedOffs is the fixed offset to
    %the experimental data
    f_out = (out.E - min(out.E))/1e9/r.h;
    mJInds = find(f_out > 3.5 & f_out < 5);
    [~,indMax] = max(out.rabi(mJInds,1));
    mjFreq_raw = f_out(mJInds(indMax));
    mJFreq = mjFreq_raw - fitOffs + fixedOffs;
end