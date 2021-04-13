clear;

Nx = 1;
Nz = 40;

waist = [623, 1064; 700, 1064]*1e-9;

t_move_ramp = [1e-3 0.1e-3];

x = 0; %linspace(0,1e-6,Nx);
z = linspace(0,3e-6,Nz);

[X,Z] = ndgrid(x,z);
xx = X(:);
zz = Z(:);

for j = 1:2
    clear gs_overlap E_out params;
    
    start_time = now;
    
    for i = 1:numel(zz)
        disp([num2str(i) '/' num2str(numel(zz))])
        separation = [xx(i) zz(i)];
        [gs_overlap(i,:),E_out(i,:),params(i)] = merge_fg_splitstep_2d_fun(separation,waist(j,:),t_move_ramp);
    end
    params = struct2table(params);
    
    save(['data_' datestr(start_time,'YYmmDD_HHMMSS') '.mat'],'params','gs_overlap','E_out')
    
    gs_overlap = reshape(permute(gs_overlap,[1 3 2]),[Nx,Nz,2]);
    
    figure(1);
    clf;
    subplot(1,2,1);
    imagesc(z*1e6,x*1e6,gs_overlap(:,:,1));
    shading flat
    xlabel('z (um)');
    xlabel('x (um)');
    colorbar
    
    subplot(1,2,2);
    imagesc(z*1e6,x*1e6,gs_overlap(:,:,2));
    shading flat
    xlabel('z (um)');
    xlabel('x (um)');
    colorbar
end