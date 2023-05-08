
if ~exist('freqsv0','var')
    [freqsv0,v0survivals,v0survivals_err] =  loadv0();
    v0offs = 691;
    freqsv0 = freqsv0 - v0offs;
end
if ~exist('freqsv22','var')
[freqsv22,v22survivals,v22survivals_err] =  loadFiles({'names_20220126_005424','names_20220126_082139',...
    'names_20220126_105135'},[1:20,44:58,64:155]);
    freqsv22 = freqsv22 - 320010.6;
    v22offs = 0;
    freqsv22 = freqsv22 - v22offs;
end
if ~exist('freqsv26','var')
    [freqsv26,v26survivals,v26survivals_err] =  loadFiles({'names_20200302_123012','names_20200302_022453','names_20200302_171940','names_20200302_175803'....
        'names_20200302_233758','names_20200302_233758','names_20200303_120922' ,'names_20200303_120922'});
    freqsv26 = freqsv26 - 325115;
    v26offs = 3;
    freqsv26 = freqsv26 - v26offs;
    [freqsv26,order] = sort(freqsv26);
    v26survivals = v26survivals(order);
    v26survivals_err = v26survivals_err(order);
end


r = constants;

r.B = (862)*1e-4; %B Field in Gauss. updated FR in model is at 864.1
r.initState = 2; %1 for FB mol, 2 for atoms
r.basis = 'aFC'; %Use fully coupled basis (F,mF,J,mJ)
r.recompute = 2;
r.Jmax = 3;%Max J to inlcude in rotational basis
r.mtot = [3:5]; %Range of mF to include in hyperfine basis
r.f_vib = 325129e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
r.pol_up = 'sigpm';
r.c3Sigma.Be = 1.139*1e9*r.h; %10.1103/PhysRevA.105.063322
r.c3Sigma.alpha1 = 0.4429*1e9*r.h; %0.340*1e9*c.h;
r.c3Sigma.alpha2 = -0.09568*1e9*r.h;
%

% figure(1)
% subplot(2,1,1)
% errorbar(freqsv0,v0survivals,v0survivals_err)
% xlim([0,16])
% subplot(2,1,2)
% plot(detuning/1e9,popn)
% xlim([0,16])
% 
% figure(1)
% errorbar(freqsv26,v26survivals,v26survivals_errs)

guessv0 = [2e-3,0,0.6,0];
guessv22 = [2e-3,3.5e9,0.6,0];
guessv26 = [2e-3,-7.25e9,0.6,0];


r.f_vib = 288690;
r.Gamma = r.h*90e6;
[residualsv0, fitparamsv0,fitcovv0, fitfunv0] = fitHf(freqsv0,v0survivals,v0survivals_err,r,guessv0);

% r.f_vib = 320000;
% r.Gamma = r.h*15e6;
% r.pol_up = 'sigpm';
% [residualsv22, fitparamsv22,fitcovv22, fitfunv22] = fitHf(freqsv22,v22survivals,v22survivals_err,r,guessv22);

r.f_vib = 325129;
r.Gamma = r.h*90e6;
r.pol_up = 'sigpm';
[residualsv26, fitparamsv26,fitcovv26, fitfunv26] = fitHf(freqsv26,v26survivals,v26survivals_err,r,guessv26);

% 
figure(1)
clf
subplot(3,1,[1:2])
fPlot = linspace(freqsv0(1),freqsv0(end),1e4);
errorbar(freqsv0,v0survivals,v0survivals_err)
hold on
plot(fPlot,fitfunv0(fitparamsv0,fPlot),'color','black');
legend({'Data','Model fit'})
xlabel('Freq [GHz]')
ylabel('Survival')
subplot(3,1,3)
scatter(freqsv0,residualsv0)
ylabel('Weighted residual')
xlabel('Freq [GHz]')

% figure(2)
% clf
% subplot(3,1,[1:2])
% fPlot = linspace(min(freqsv22),max(freqsv22),1e4);
% errorbar(freqsv22,v22survivals,v22survivals_err)
% hold on
% plot(fPlot,fitfunv22(fitparamsv22,fPlot),'color','black');
% legend({'Data','Model fit'})
% xlabel('Freq [GHz]')
% ylabel('Survival')
% subplot(3,1,3)
% scatter(freqsv22,residualsv22)
% ylabel('Weighted residual')
% xlabel('Freq [GHz]')
% 
figure(3)
clf
subplot(3,1,[1:2])
fPlot = linspace(min(freqsv26),max(freqsv26),1e4);
errorbar(freqsv26,v26survivals,v26survivals_err)
hold on
plot(fPlot,fitfunv26(fitparamsv26,fPlot),'color','black');
legend({'Data','Model fit'})
xlabel('Freq [GHz]')
ylabel('Survival')
subplot(3,1,3)
scatter(freqsv26,residualsv26)
ylabel('Weighted residual')
xlabel('Freq [GHz]')

function [residuals, fitparams,fitcov, fitfun] = fitHf(freqs,survivals,survivals_err,ramanSet,guess)
    % %Compute effective Hamiltonian for up-leg (and down-leg) including
    % %rovibrational contribution to TDM
    raman_data = raman_effective_hamiltonian(ramanSet.B,ramanSet.basis,ramanSet.recompute,ramanSet.Jmax,ramanSet.mtot,ramanSet.f_vib);

    fitfun = @(b,x) findRsc(x,raman_data,ramanSet,ramanSet.initState,ramanSet.pol_up,b(1),b(2),b(3),b(4),ramanSet.c3Sigma.Gamma);
    [fitparams,residuals,~,fitcov] = nlinfit(freqs, survivals,fitfun,guess,'weights',1./survivals_err.^2);
end

function popn = findRsc(f,raman_data,const,initState,pol_up,t,f_offs,scale,y_offs,Gamma)
    % initial state
    psi_init = zeros(size(raman_data.f.E));
    %Select lowest energy state fromfeshbach hamiltonian, corresponding to the
    %FB molecule
    if initState == 1
        indState = find(raman_data.f.E < raman_data.f.E_lowest_chan_threshold);
        psi_init(indState(end)) = 1; %Find highest energy state below threshold
    elseif initState == 2
        indState = find(raman_data.f.E > raman_data.f.E_lowest_chan_threshold);
        psi_init(indState(1)) = 1; %Find lowest energy state above threshold
    end    

% laser stuff for figuring out realistic scattering rates
    power = 100e-3; %W
    waist = 50e-6; %m
    p_up = sphten(polSwitch(pol_up));

    Efield = sqrt(4*const.eta0*power/(pi*waist^2));
    freq_offs = min(min(raman_data.c.E - raman_data.f.E))/const.h + f_offs;

    %% 
    %Calculate up-leg rabi rates to each state and corresponding total
    %scattering rates
%     Gamma = const.h*90e6;
    rabi = Efield*sum(raman_data.H_up.*psi_init.*permute(p_up,[1 3 4 2]),4);
    delta = raman_data.c.E - raman_data.f.E - const.h*reshape(freq_offs + f*1e9,1,1,[]);
    Rsc = permute(sum(sum((Gamma/(2*const.h)) * 2*abs(rabi).^2./(2*abs(rabi).^2 + 4*delta.^2 + Gamma.^2),2),1),[2 3 1]);
    popn = abs(scale)*exp(-abs(Rsc*t)) + y_offs;
%     popn = exp(-Rsc*t);
end

function [freqsv0,v0survivals,v0survivals_err] =  loadv0()
    paScan = {'20180417_234100','20180418_042329'};
    datadir = 'C:\NaCs1pt5_Data_temp\Data\';
    loadPA{1} =  load([datadir, paScan{1}(1:8),'\', 'VPAScan_', paScan{1},'.mat']);
    nameFilesv0{1} = loadPA{1}.filelist;
    loadPA{2} =  load([datadir, paScan{2}(1:8),'\', 'VPAScan_', paScan{2},'.mat']);
    nameFilesv0{2} = loadPA{2}.filelist;
    nameFilesv0 = cat(2,nameFilesv0{1},nameFilesv0{2});
    freqsv0{1} = loadPA{1}.fWavemeter;
    freqsv0{2} = loadPA{2}.fWavemeter;
    freqsv0 = cat(2,freqsv0{1},freqsv0{2});
    freqsv0 = freqsv0(1:334);
    v0survivals = [];
    v0survivals_err = [];
    for i = 1:334
        loaddat =  load([datadir, nameFilesv0{i}(1:8),'\', 'data_', nameFilesv0{i},'.mat']);
        if i == length(nameFilesv0)
            disp(1)
        end
        singleAtomLogical = loaddat.Analysis.SingleAtomLogical;
        loadLogical = sum(squeeze(singleAtomLogical(1:2,:)),1) == 2;
        survLogical = sum(squeeze(singleAtomLogical(3:4,:)),1) == 2;
        [v0survivals(i), v0survivals_err(i)] = find_survival(survLogical, loadLogical,loadLogical, [1], 1);
    end

    freqs{1} = freqsv0(1:end-1) + 288000;
    survivals{1} = v0survivals;
    survivals_errs{1} = v0survivals_err;
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