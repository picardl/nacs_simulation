clear;

% files = {'data_212112_220452.mat','data_210612_230455.mat'};
% files = {'data_210414_114847.mat'};
% files = {'data_210414_114847.mat','data_210414_115805.mat'};
% files = {'data_210414_133931.mat','data_210414_141940.mat'}; % phi = 5 deg
% % files = {'data_210414_133947.mat','data_210414_142017.mat'}; % phi = 0 deg

% files = {'data_210414_133931.mat','data_210414_142017.mat'}; % phi = 0 deg
files = {'data_210414_133947.mat','data_210414_141940.mat'}; % phi = 5 deg

% files = {'data_210414_161602.mat','data_210414_162923.mat'}; % 5 deg
% files = {'data_210414_161619.mat','data_210414_162930.mat'}; % 0 deg

N = 20;

for i = 1:numel(files)
    data = load(files{i});
%     z = data.z;
%     x = data.x;
    z = reshape(data.params.align_err(:,1),N,N);
    x = reshape(data.params.align_err(:,2),N,N);
    waist = mean(data.params.waist,1)*1e9;
    
    gs_overlap = reshape(permute(data.gs_overlap,[1 3 2]),N,N,2);
    E_out = reshape(permute(data.E_out,[1 3 2]),N,N,2);
    
    figure(100*i+1);
    subplot(1,2,1);
    imagesc(x(1,:)'*1e6,z(:,1)*1e6,gs_overlap(:,:,1));
    colorbar;
    title(['Na g.s. overlap, waists = ' num2str(waist(1)) ',' num2str(waist(2)) ' nm'])
    xlabel('axial (um)')
    ylabel('radial (um)')
    caxis([0 1]);
    subplot(1,2,2);
    imagesc(x(1,:)'*1e6,z(:,1)*1e6,gs_overlap(:,:,2));
    colorbar;
    title(['Cs g.s. overlap, waists = ' num2str(waist(1)) ',' num2str(waist(2)) ' nm'])
    xlabel('axial (um)')
    ylabel('radial (um)')
    set(gcf,'units','inches');
    p = get(gcf,'position');
    set(gcf,'units','inches','position',[p(1) p(2) 8 3])
    caxis([0 1]);
    
    figure(100*i+2);
    subplot(1,2,1);
    imagesc(x(1,:)'*1e6,z(:,1)*1e6,(E_out(:,:,1)-1)/2);
    colorbar;
     title(['Na nbar, waists = ' num2str(waist(1)) ',' num2str(waist(2)) ' nm'])
     xlabel('axial (um)')
    ylabel('radial (um)')
%     caxis([0 1]);
    subplot(1,2,2);
    imagesc(x(1,:)'*1e6,z(:,1)*1e6,(E_out(:,:,2)-1)/2);
    colorbar;
   title(['Cs nbar, waists = ' num2str(waist(1)) ',' num2str(waist(2)) ' nm'])
    xlabel('axial (um)')
    ylabel('radial (um)')
    set(gcf,'units','inches');
    p = get(gcf,'position');
    set(gcf,'units','inches','position',[p(1) p(2) 8 3])
%     caxis([0 1]);
end