clear;

% files = {'data_212112_220452.mat','data_210612_230455.mat'};
% files = {'data_210413_101955.mat','data_210413_102017.mat'};
% files = {'data_210413_105841.mat','data_210413_110007.mat'};
files = {'data_210413_114936.mat','data_210413_115451.mat'};


for i = 1:2
    data = load(files{i});
%     z = data.z;
%     x = data.x;
%     x = reshape(data.params.align_err(:,1),10,10);
%     z = reshape(data.params.align_err(:,2),10,10);
    z = data.params.align_err(:,2);
    waist = mean(data.params.waist(:,1)*1e9);
    gs_overlap = permute(data.gs_overlap,[1 3 2]);
    E_out = permute(data.E_out,[1 3 2]);
    
    figure(1+100*i);
%     clf;
    hold on; box on;
    for j = 1:2
        plot(z*1e6,gs_overlap(:,:,j))
    end
    hold off;
    xlabel('axial (um)')
    ylabel('g.s. wfn overlap')
    legend(['na, waist = ' num2str(waist)],['cs, waist = ' num2str(waist)])
%     subplot(1,2,1);
%     pcolor(z,x,gs_overlap(:,:,1));
%     shading flat
%     colorbar;
%     subplot(1,2,2);
%     pcolor(z,x,gs_overlap(:,:,2));
%     shading flat
%     colorbar;
    

    figure(2+100*i);
%     clf;
    hold on; box on;
    for j = 1:2
        plot(z*1e6,E_out(:,:,j))
    end
    hold off;
    xlabel('axial (um)')
    ylabel('E/(\omega/2)')
    legend(['na, waist = ' num2str(waist)],['cs, waist = ' num2str(waist)])
%     subplot(1,2,1);
%     pcolor(z,x,E_out(:,:,1));
%     shading flat
%     colorbar;
%     subplot(1,2,2);
%     pcolor(z,x,E_out(:,:,2));
%     shading flat
%     colorbar;
end