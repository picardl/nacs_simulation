
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

r.B = (864)*1e-4; %B Field in Gauss. updated FR in model is at 864.1
r.initState = 2; %1 for FB mol, 2 for atoms
r.basis = 'aFC';
r.recompute = 2;
r.Jmax = 3;%Max J to inlcude in rotational basis
r.mtot = [2:6]; %Range of mF to include in hyperfine basis
r.f_vib = 320000e9; %Freq of transition in Hz. Nearest vibrational level to this value will be selected.
r.pol_up = 'sigpm';
r.c3Sigma.Be = 1.04*1e9*r.h; %10.1103/PhysRevA.105.063322
r.c3Sigma.alpha1 =  0.415*1e9*r.h; %0.340*1e9*c.h;
r.c3Sigma.alpha2 =  0.1125*1e9*r.h;
r.c3Sigma.Gamma = 15e6*r.h;
r.power = 100e-3;
r.waist = 50e-6;
r.N_tot = 0;
vib_data = load('../data/cbB_210705_125650.mat');
E_vib = vib_data.out.E*r.hartree + r.Cs_D12_weighted;
ind = find(abs(E_vib-r.f_vib*r.h) == min(abs(E_vib-r.f_vib*r.h)));
r.E_vib = E_vib(ind);

f_off = 3;
pre = precomp_ops(r);
% pre = precomp_hp(r);
% pre.psi_init = zeros(size(pre.psi_init));
% pre.psi_init(4) = 1;

guessv0 = [2e-3,0,0.6,0];
guessv22 = [r.c3Sigma.Be/(1e9*r.h),r.c3Sigma.alpha1/(1e9*r.h),r.c3Sigma.alpha2/(1e9*r.h),f_off,1e-3,0.55,0.075];
guessv26 = [2e-3,-7.25e9,0.6,0.01];

fitfun = @(b,x) HamnRSC(x,r,pre,b(1),b(2),b(3),b(4),guessv22(5),guessv22(6),guessv22(7));
fitfun_t = @(b,x) HamnRSC(x,r,pre,b(1),b(2),b(3),b(4),b(5),guessv22(6),guessv22(7));
lsqfitfun =@(b) sum((v22survivals - HamnRSC(freqsv22,r,pre,b(1),b(2),b(3),b(4),guessv22(5),guessv22(6),guessv22(7))).^2);
lsqfitfun_t = @(b) sum((v22survivals - HamnRSC(freqsv22,r,pre,b(1),b(2),b(3),b(4),b(5),guessv22(6),guessv22(7))).^2);

NOff = 1;
off_range = linspace(2.7,2.9,NOff);

BeScan = [0.91];
alpha1Scan = [0.3:0.1:0.51];
alpha2Scan = [0.07:0.05:0.12];

ScanPts = length(BeScan)*length(alpha1Scan)*length(alpha2Scan);

% resids = zeros(length(BeScan),length(alpha1Scan),length(alpha2Scan));
% minOffs = zeros(length(BeScan),length(alpha1Scan),length(alpha2Scan));
% f = waitbar(0,"Scan progress");
% BeLin = zeros(1,ScanPts);
% a1Lin = zeros(1,ScanPts);
% a2Lin = zeros(1,ScanPts);
% residsLin = zeros(1,ScanPts);
% minOffsLin = zeros(1,ScanPts);
% indLin = 1;

% disp(['Estimated scan time: ',num2str(ScanPts*NOff*0.29/60),' min'])
% 
% try
%     for i = 1:length(BeScan)
%         for j = 1:length(alpha1Scan)
%             for k = 1:length(alpha2Scan)
%                 residsOff = zeros(1,NOff);
%                 for l = 1:NOff
%                     residsOff(l) = lsqfitfun([BeScan(i),alpha1Scan(j),alpha2Scan(k),off_range(l)]);
%                 end
%               [moff,indmin] = min(residsOff);
%               resids(i,j,k) = moff;
%               minOffs(i,j,k) = off_range(indmin);
%               BeLin(indLin) = BeScan(i);
%               a1Lin(indLin) = alpha1Scan(j);
%               a2Lin(indLin) = alpha2Scan(k);
%               residsLin(indLin) = resids(i,j,k);
%               minOffsLin(indLin) = minOffs(i,j,k);
%               waitbar((i-1)/length(BeScan) + (j-1)/length(alpha1Scan)/length(BeScan) + k/length(alpha2Scan)/length(BeScan)/length(alpha1Scan),f);
%               indLin = indLin + 1;
%             end
%          end
%     end
%     close(f);
% catch
%     disp('error in scan, continuing to plot');
% end

% saveDat = struct();
% saveDat.resids = resids;
% saveDat.minOffs = minOffs;
% saveDat.BeLin = BeLin;
% saveDat.a1Lin = a1Lin;
% saveDat.a2Lin = a2Lin;
% saveDat.residsLin= residsLin;
% saveDat.minOffsLin =minOffsLin;
% saveDat.r = r;
% fn = ['../data/c3Sigma_hfFitData_',datestr(now,'YYmmDD_HHMMSS'), '.mat'];
% save(fn,'saveDat')
% disp(fn);

% [X,Y,Z] = meshgrid(BeScan,alpha1Scan,alpha2Scan);
% slice(X,Y,Z, resids,[1],[0.4],[0.1])    % display the slices
% cb = colorbar;                                  % create and label the colorbar
% cb.Label.String = 'Residual';

% figure(2)
% scatter3(BeLin,a1Lin,a2Lin,200,residsLin,'filled')    % draw the scatter plot
% ax = gca;
% ax.XDir = 'reverse';
% view(-31,14)
% xlabel('Be')
% ylabel('Alpha_{Na}')
% zlabel('Alpha_{Cs}')
% cb = colorbar;                                     % create and label the colorbar

figure(3)
clf
plot(freqsv22,v22survivals)
hold on
plot(freqsv22,fitfun_t([0.9893    0.2280    0.1443    2.8222    0.00104],freqsv22))
lsqfitfun_t([0.9893    0.2280    0.1443    2.8222    0.00104])
xlabel('Freq [GHz]')
ylabel('Population')

% [X,Y,Z] = meshgrid(BeScan,alpha1Scan,alpha2Scan);
% slice(X,Y,Z, resids,[1],[0.4],[0.1])    % display the slices
% cb = colorbar;                                  % create and label the colorbar
% cb.Label.String = 'Residual';

% figure(6)
% clf
% [X,Y] = meshgrid(alpha1Scan,alpha2Scan);
% surf(X,Y,squeeze(resids(1,:,:))')
% xlabel('Alpha_{Na}')
% ylabel('Alpha_{Cs}')
% c = colorbar();
% c.Label.String = 'Residuals';
% zlabel('Residuals')

if 0
    %%
    options = optimset;
    opstat = statset;
    opstat.Display = 'iter';
    options.Display = 'iter';
    fitfun2 = @(b,x) fitfun_t([b(1),0.33,0.0825,b(2),0.55e-3],x);
    lsqfitfun2 = @(b) lsqfitfun_t([b(1),b(2),b(3),b(4),b(5)]);
%     [beta,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(freqsv22,v22survivals,fitfun2,[0.904,2.8],opstat,'weights',1./v22survivals_err.^2);

%     BRange = linspace(0.89,0.93,100);
%     lsq = [];
%     for i = 1:100
%         lsq(i) = lsqfitfun2([BRange(i),2.8053]);
%     end
%     figure(12)
%     plot(BRange,lsq)
    fitvals = fminsearch(lsqfitfun2,[0.97    0.3564    0.0972    2.7989    0.003],options)
end
if 0
   %% unpack save data
    resids = saveDat.resids;
    minOffs = saveDat.minOffs;
    BeLin = saveDat.BeLin;
    a1Lin = saveDat.a1Lin;
    a2Lin = saveDat.a2Lin;
    residsLin = saveDat.residsLin;
    minOffsLin = saveDat.minOffsLin;
    r = saveDat.r;
end

function  popn = HamnRSC(f,r,pre,Be,alpha1,alpha2,f_offs,t,scale,y_offs)
    % %Compute effective Hamiltonian for up-leg (and down-leg) including
    % %rovibrational contribution to TDM

    r.c3Sigma.Be = Be*1e9*r.h;
    r.c3Sigma.alpha1 = alpha1*1e9*r.h;
    r.c3Sigma.alpha2 = alpha2*1e9*r.h;
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
    popn = abs(scale)*exp(-abs(Rsc*t)) + min(y_offs,0.15);
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