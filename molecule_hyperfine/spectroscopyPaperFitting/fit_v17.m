[freqs,survivals,survivals_errs] =  loadFiles({'names_20220112_190342','names_20220113_011604'},[1:2:232]);
offs = 313346.3;
freqs = freqs - offs + 80e-3;

r = constants;
% load('../data/precomp_all_v.mat');
pre = prev17;

r.B = 860*1e-4; %B Field in Gauss. updated FR in model is at 864.1
r.initState = 2; %1 for FB mol, 2 for atoms
r.basis = 'aFC';
r.recompute = 2;
r.Jmax = 3;%Max J to inlcude in rotational basis
r.mtot = [3:6]; %Range of mF to include in hyperfine basis
r.f_vib = 313346.3e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
r.pol_up = 'sigpm';
r.pol_frac = 0.5;
r.c3Sigma.Be = 1.04*1e9*r.h; %10.1103/PhysRevA.105.063322
r.c3Sigma.alpha1 =  0.7*1e9*r.h; %0.340*1e9*c.h;
r.c3Sigma.alpha2 =  0.09*1e9*r.h;
r.c3Sigma.Gamma = 100e6*r.h;
r.power = 100e-3;
r.waist = 30e-6;
r.N_tot = 0;
vib_data = load('../data/cbB_210703_014725.mat');
E_vib = vib_data.out.E*r.hartree + r.Cs_D12_weighted;
ind = find(abs(E_vib-r.f_vib*r.h) == min(abs(E_vib-r.f_vib*r.h)));
r.E_vib = E_vib(ind);

f_off = 3;
% pre = precomp_ops(r);
pre.psi_init = zeros(size(pre.psi_init));
pre.psi_init(4) = 1;

guess = [r.c3Sigma.Be/(1e9*r.h),r.c3Sigma.alpha1/(1e9*r.h),r.c3Sigma.alpha2/(1e9*r.h),f_off,1e-3,0.5,mean(survivals(1:10))-min(survivals),min(survivals)];
fitfun = @(b,x) HamnRSC(x,r,pre,b(1),b(2),b(3),b(4),b(5),guess(6),guess(7),guess(8));
fitfunNoCs = @(b,x) HamnRSC(x,r,pre,b(1),b(2),0.05028827,b(3),b(4),guess(6),guess(7),guess(8));
lsqfitfun =@(b) sum((survivals - HamnRSC(freqs,r,pre,b(1),b(2),b(3),b(4),b(5),guess(6),guess(7),guess(8))).^2);
lsqfitfun2 =@(b) sum((survivals - HamnRSC(freqs,r,pre,b(1),b(2),guess(3),b(3),guess(5),guess(6),guess(7),guess(8))).^2);
lsqfitfunNoCs =@(b) sum((survivals - HamnRSC(freqs,r,pre,b(1),b(2),0.050288,b(3),b(4),guess(6),guess(7),guess(8))).^2);

% 
% % 
figure(3)
clf
ftfreqs = linspace(min(freqs),max(freqs),1000);
plot(freqs,survivals)
hold on
plot(ftfreqs,fitfunNoCs([1.0036    0.5741   -0.2160    0.0850],ftfreqs))
lsqfitfunNoCs([1.0036    0.5741   -0.2160    0.0850])
xlabel('Freq [GHz]')
ylabel('Population')

if 0
    %%
    r.c3Sigma.Be =1.0036*1e9*r.h;
    r.c3Sigma.alpha1 = 0.5741*1e9*r.h;
    r.c3Sigma.alpha2 =  0.05028827*1e9*r.h;
    r.pol_up = 0.9;
    out = upleg_hamn_precomp(pre,r,r.B); 
    mJFreq = offToMj1(r,out,-0.2160,offs)
end

if 0
    %%
    options = optimset;
    opstat = statset;
    opstat.Display = 'iter';
    options.Display = 'iter';
%     fitfun2 = @(b,x) fitfun_t([b(1),0.33,0.0825,b(2),0.55e-3],x);
%     lsqfitfun2 = @(b) lsqfitfun_t([b(1),b(2),b(3),b(4),b(5)]);

%     fitvals = fminsearch(lsqfitfun,[1   0.7105   -0.0253   0.6    0.5],options)
%     fitvals = fminsearch(lsqfitfunNoCs,[1   0.7105   0.6    0.5],options)

    [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(freqs,survivals,fitfunNoCs,[1.0091    0.5779   -0.1895    0.0830],opstat,'weights',1./survivals_errs.^2);

end

if 0
   %% 
    a1 = [0.3:0.02:0.8];
    a2 = [0.01:0.005:0.1];
    lsq = [];
    for i = 1:length(a1)
        for j = 1:length(a2)
           lsq(i,j) =  lsqfitfun([1.0036    a1(i)   a2(j)   -0.2160    0.0850])
        end
    end
    [X,Y] = ndgrid(a1,a2);
    figure(11)
    surf(X,Y,lsq)
    xlabel('alpha Na [GHz]')
    ylabel('alpha Cs [GHz]')
end

function  popn = HamnRSC(f,r,pre,Be,alpha1,alpha2,f_offs,t,pol_up,scale,y_offs)
    % %Compute effective Hamiltonian for up-leg (and down-leg) including
    % %rovibrational contribution to TDM

    r.c3Sigma.Be = Be*1e9*r.h;
    r.c3Sigma.alpha1 = alpha1*1e9*r.h;
    r.c3Sigma.alpha2 = alpha2*1e9*r.h;
    r.pol_up = pol_up;
    out = upleg_hamn_precomp(pre,r,r.B);
    popn = findRsc(f,out.E,pre.f.E,out.rabi,r,t,f_offs,scale,y_offs,r.c3Sigma.Gamma);
end

function popn = findRsc(f,E,E_fb,rabi,const,t,f_offs,scale,y_offs,Gamma)
    %% 
    %Calculate up-leg rabi rates to each state and corresponding total
    %scattering rates
    freq_offs = min(min(E - E_fb))/const.h + f_offs*1e9;
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
            survivals(:,i,:) = permute(loaddat.Analysis.SurvivalProbability,[1 3 2]);
            survivals_errs(:,i,:) = permute(loaddat.Analysis.SurvivalProbabilityErr,[1 3 2]);
        catch
            try
            survivals(1:5,i,:) = permute(loaddat.Analysis.SurvivalProbability,[1 3 2]);
            survivals_errs(1:5,i,:) = permute(loaddat.Analysis.SurvivalProbabilityErr,[1 3 2]);
            catch
                try
                    survivals(1:5,i,:) = permute(loaddat.Analysis.SurvivalProbability(:,2),[1 3 2]);
                    survivals_errs(1:5,i,:) = permute(loaddat.Analysis.SurvivalProbabilityErr(:,2),[1 3 2]);
                catch
                    try
                        survivals(1:5,i,:) = loaddat.Analysis.SurvivalProbability([3,1,2,4,5],end);
                        survivals_errs(1:5,i,:) =loaddat.Analysis.SurvivalProbabilityErr([3,1,2,4,5],end);
                    catch
                        disp(['error loading file ' nameFiles{i}]);
                        survivals(:,i,:) = nan;
                        survivals_errs(:,i,:) = nan;
                    end
                end
            end   
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
    survivals = survivals(1,order,:);
    survivals_errs = survivals_errs(1,order,:);
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
    mJInds = find(f_out > 2.5 & f_out < 4.5);
    [~,indMax] = max(out.rabi(mJInds,1));
    mjFreq_raw = f_out(mJInds(indMax));
    mJFreq = mjFreq_raw - fitOffs + fixedOffs;
end