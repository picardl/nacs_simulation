clear;

files = {'data_210414_113954.mat','data_210414_114208.mat'};
files = {'data_210414_114632.mat','data_210414_114618.mat'};
files = {'data_210414_114957.mat','data_210414_114920.mat'};

for i = 1:numel(files)
    data = load(files{i});
    pow = data.params.power(:,1);
    waist = mean(data.params.waist(:,1)*1e9);
    gs_overlap = permute(data.gs_overlap,[1 3 2]);
    E_out = permute(data.E_out,[1 3 2]);
    
    figure(1+100*i);
    hold on; box on;
    for j = 1:2
        plot(pow,gs_overlap(:,:,j))
    end
    hold off;
    xlabel('Na power (mW)')
    ylabel('g.s. wfn overlap')
    legend(['na, waist = ' num2str(waist)],['cs, waist = ' num2str(waist)])
    

    figure(2+100*i);
    hold on; box on;
    for j = 1:2
        plot(pow,E_out(:,:,j))
    end
    hold off;
    xlabel('Na power (mW)')
    ylabel('E/(\omega/2)')
    legend(['na, waist = ' num2str(waist)],['cs, waist = ' num2str(waist)])
end