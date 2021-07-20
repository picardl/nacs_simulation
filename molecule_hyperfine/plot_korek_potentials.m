function plot_korek_potentials

const = constants();

%% work out energy offset
% Eoffs_guess = -0.331957;
%
% new_data = readtable(['../lib/korek/NaCs1'],'NumHeaderLines',1);
%
% R = new_data{:,1};
% V = new_data{:,2};
%
% C6 = 1.8353;
% C8 = 3.7927;
%
% Vlr = - C6.*(R*const.abohr*1e10).^-6 - C8.*(R*const.abohr*1e10).^-8;
%
% plot(R,V,R,Vlr+Eoffs_guess);

%% re-save data in individual files
% states_count = 0;
%
% for i = 1:6
%     new_data = readtable(['../lib/korek/NaCs' num2str(i)],'NumHeaderLines',1);
%
%     colnames = cellfun(@(x) ['var' num2str(x)],num2cell((states_count+1):(states_count+size(new_data,2)-1)),'un',false);
%     colnames = ['R' colnames];
%     new_data.Properties.VariableNames = colnames;
%
%     new_data = unique(new_data,'rows');
%
%     for j = 2:size(new_data,2)
%         save_data = table;
%         save_data.R = new_data.R;
%         save_data.V = new_data{:,j};
%
% %         writetable(save_data,['../lib/korek/state' num2str(states_count+j-1) '.xlsx']);
%     end
%
%     states_count = states_count + size(new_data,2)-1;
%
% end

% B_R = const.hbar^2./(2*const.mu_nacs*(R*const.abohr).^2);
% hold on
% plot(R,B_R/const.hartree)
% hold off

% figure(1);
% clf;
% so_files = {'SO0-','SO0+','SO2','SO3','SO1'};
% for i = 1:numel(so_files)
%     data = readtable(['../lib/korek/' so_files{i}],'NumHeaderLines',1);
%     hold on;
%     for j = 2:size(data,2)
%         plot(data{:,1},data{:,j},'-')
%     end
%     hold off;
% end
%
% r = data{:,1};
% V = data{:,2};
%
% rfit = 5.3;
%
% Vfun =@(b,r) 1./r + b(1)./r.^4 + b(2);
%
% b = nlinfit(r(r<rfit),V(r<rfit),Vfun,[1 1]);
%
% hold on;
% plot(r,V,'-k',r,Vfun(b,r),'--k');
% hold off;
%
% % ylim([-0.36 -0.26])
%
% b

%%

% col = distinguishable_colors(6);
% 
% R = linspace(4.5,20,1e2);
% 
% terms = {'1S','3S','1P','3P','1D','3D'};
% 
% figure(1);
% clf;
% for i = 1:numel(terms)
%     j = 1;
%     while true
%         try
%             V = korek_potential(R,terms(i),j);
%             hold on; box on;
%             p = plot(R,V,'color',col(i,:));
%             hold off;
%             if j==1
%                 pp(i) = p;
%                 leg(i) = terms(i);
%             end
%             j = j+1;
%         catch
%             break
%         end
%     end
%     %     data = readtable(['../lib/korek/state' num2str(i) '.xlsx']);
%     
% end
% legend(pp,leg);


%%
Rtest = 7.2;
terms = {'1S','3S','1P','3P','1D','3D'};
for i = 1:6
    new_data = readtable(['../lib/korek/NaCs' terms{i} '.csv'],'headerlines',1);
    cols = size(new_data,2);
    col_names = strcat(terms(i),'_',arrayfun(@num2str,1:(cols-1),'un',0));
    new_data.Properties.VariableNames = ['R',col_names];
    
    
    if i==1
        data = new_data;
    else
        for j = 2:size(new_data,2)
            x = new_data{:,1};
            y = new_data{:,j};
            if ~(numel(unique(x))==numel(x))
                ux = unique(x);
                [~,subs] = ismember(x,ux);
                uy = accumarray(subs,y,[],@mean);
                x = ux;
                y = uy;
            end
            cut = isnan(y);
            x(cut) = [];
            y(cut) = [];
            [x,order] = sort(x);
            y = y(order);
            data{:,col_names{j-1}} = interp1(x,y,data.R);
        end
%         data = outerjoin(data,new_data,'mergekeys',true);
    end
end
[~,order] = sort(data{data.R==Rtest,2:end});
data = data(:,[1 order+1]);

alphX = 'XABCDEFGHIJKLMNOPQRSTUVWYZ';
alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';

colnames = data.Properties.VariableNames(2:end);
colnames = cellfun(@(x) x(1:end-2),colnames,'un',0);

singlet = cellfun(@(x) str2double(x(1))==1,data.Properties.VariableNames(2:end));
triplet = cellfun(@(x) str2double(x(1))==3,data.Properties.VariableNames(2:end));

singlet_labels = cellstr(alphX(1:numel(colnames(singlet)))');
triplet_labels = cellstr(lower(alph(1:numel(colnames(triplet))))');

colnames(singlet) = strcat(singlet_labels',colnames(singlet));
colnames(triplet) = strcat(triplet_labels',colnames(triplet));

data.Properties.VariableNames = ['R',colnames];

% figure(1);
% clf
% hold on
% for i = 2:size(data,2)
%     plot(data.R,data{:,i})
% end
% hold off
% legend(data.Properties.VariableNames(2:end))

writetable(data,'../lib/korek/korek_data.xlsx');


end