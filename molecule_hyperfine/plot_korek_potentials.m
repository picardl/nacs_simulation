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

col = distinguishable_colors(6);

R = linspace(4.5,20,1e2);

terms = {'1S','3S','1P','3P','1D','3D'};

figure(1);
clf;
for i = 1:numel(terms)
    j = 1;
    while true
        try
            V = korek_potential(R,terms(i),j);
            hold on; box on;
            p = plot(R,V,'color',col(i,:));
            hold off;
            if j==1
                pp(i) = p;
                leg(i) = terms(i);
            end
            j = j+1;
        catch
            break
        end
    end
    %     data = readtable(['../lib/korek/state' num2str(i) '.xlsx']);
    
end
legend(pp,leg);





end