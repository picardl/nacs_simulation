clear;

files = {'data_212112_220452.mat','data_210612_230455.mat'};

for i = 1:2
    data = load(files{i});
    z = data.z;
    x = data.x;
    
    gs_overlap = reshape(permute(data.gs_overlap,[1 3 2]),numel(z),numel(x),2);
    E_out = reshape(permute(data.E_out,[1 3 2]),numel(z),numel(x),2);
    
    figure(100*i+1);
    subplot(1,2,1);
    imagesc(z*1e6,x*1e6,gs_overlap(:,:,1));
    colorbar;
    title('Na ground state overlap')
    xlabel('axial (um)')
    ylabel('radial (um)')
    subplot(1,2,2);
    imagesc(z*1e6,x*1e6,gs_overlap(:,:,2));
    colorbar;
    title('Cs ground state overlap');
    xlabel('axial (um)')
    ylabel('radial (um)')
    
    figure(100*i+2);
    subplot(1,2,1);
    imagesc(z*1e6,x*1e6,log10(E_out(:,:,1)));
    colorbar;
     title('Na log10(energy)')
     xlabel('axial (um)')
    ylabel('radial (um)')
    subplot(1,2,2);
    imagesc(z*1e6,x*1e6,log10(E_out(:,:,2)));
    colorbar;
    title('Cs log10(energy)');
    xlabel('axial (um)')
    ylabel('radial (um)')
end