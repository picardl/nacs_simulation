function plot_korek_potentials

const = constants();

%% work out energy offset
Eoffs_guess = -0.331957;

new_data = readtable(['../lib/korek/NaCs1'],'NumHeaderLines',1);

R = new_data{:,1};
V = new_data{:,2};

C6 = 1.8353;
C8 = 3.7927;

Vlr = (C6/R).^6 - (C8./R).^8;

plot(R,V,R,Vlr+Eoffs_guess);

%% re-save data in individual files
states_count = 0;

for i = 1:6
    new_data = readtable(['../lib/korek/NaCs' num2str(i)],'NumHeaderLines',1);
    
    colnames = cellfun(@(x) ['var' num2str(x)],num2cell((states_count+1):(states_count+size(new_data,2)-1)),'un',false);
    colnames = ['R' colnames];
    new_data.Properties.VariableNames = colnames;
    
    new_data = unique(new_data,'rows');
    
    for j = 2:size(new_data,2)
        save_data = table;
        save_data.R = new_data.R;
        save_data.V = new_data{:,j};
        
        writetable(save_data,['../lib/korek/state' num2str(states_count+j-1) '.xlsx']);
    end
    
    states_count = states_count + size(new_data,2)-1;
    
end

% B_R = const.hbar^2./(2*const.mu_nacs*(R*const.abohr).^2);
% hold on
% plot(R,B_R/const.hartree)
% hold off


% so_files = {'SO0-','SO0+','SO1','SO2','SO3'};
% for i = 1:numel(so_files)
%     data = readtable(['../lib/korek/' so_files{i}],'NumHeaderLines',1);
%     hold on;
%     for j = 2:size(data,2)
%         plot(data{:,1},data{:,j},'-')
%     end
%     hold off;
% end

end