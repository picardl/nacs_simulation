clear;

data=load('../data/d_all.mat');

B = [0 600 433 860];

figure(1);
clf;
hold on;
for i = 1:numel(data.freqs_all)
    errorbar(data.freqs_all{i},data.survivals_all{i}(3,:)+B(i)/300,data.survivals_errs_all{i}(3,:),'.-');
end
hold off;
xlim([325125 325135])

% y = [];
% dy = [];
% for i = 1:numel(names)
%     
%     data = load(['../data/20200802/data_' names{i} '.mat']);
%     
%     y = cat(1,y,data.Analysis.SurvivalProbability);
%     dy = cat(1,dy,data.Analysis.SurvivalProbabilityErr);
%     
% end
% 
% y = reshape(y,7,[]);
% 
% plot(y(3,:)')