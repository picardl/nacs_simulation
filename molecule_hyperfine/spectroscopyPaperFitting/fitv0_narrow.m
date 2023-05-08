 [freqs,survivals,survivals_errs] = loadFiles('20190503_pa_names2');
 freqs = freqs - 306000;

 

fitfun =@(b,x) b(1)*(0.5*b(3))^2./((0.5*b(3)).^2+(x-b(2)).^2) + b(4);
guess = [-0.4 698.6 0.05 0.3];
[beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(freqs,survivals,fitfun,guess,'weights',1./survivals_errs.^2);

f_ovsmpl = linspace(min(freqs),max(freqs),1000);

 figure(1)
 clf
 errorbar(freqs,survivals,survivals_errs)
 hold on
 plot(f_ovsmpl,fitfun(beta,f_ovsmpl))
 xlabel('Frequency 288xxx[GHz]')
 ylabel('Na + Cs survival')

function [freqsv0,v0survivals,v0survivals_err] =  loadv0()
    datadir = 'C:\Users\lewis\OneDrive\Documents\Current Work\Ni Lab\Papers\NaCsSpectroscopy\Data\20190503';
    loadPA =  load([datadir, '\data_20190503_pa_names2.mat']);
    nameFilesv0 = loadPA.names;
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
    datadir = 'C:\Users\lewis\OneDrive\Documents\Current Work\Ni Lab\Papers\NaCsSpectroscopy\Data\';    
    
    
    if ~iscell(namesfile)
        try
            nameFiles = load([datadir namesfile(1:8) '\' namesfile '.mat']).names;
        catch
            nameFiles =  load([datadir, namesfile(1:8),'\', 'data_', namesfile,'.mat']); % load(FindDataFile(nameFiles{i})); % load([datadir{j},'data_', nameFiles{i},'.mat']);
        end
    else
        nameFiles = {};
        for i = 1:length(namesfile)
            try
                nameFiles = [nameFiles, load([datadir namesfile{i}(1:8) '\' namesfile{i} '.mat']).names];
            catch
                nameFiles = [nameFiles, load([datadir namesfile{i}(1:8) '\' namesfile{i} '.mat']).names];
            end
        end
    end 
    if ~exist('inds','var')
        inds = 1:length(nameFiles);
    end
    nameFiles = nameFiles.names;
%     nameFiles = nameFiles(inds);
    
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