clear;

Nx = 11;
Nz = 1;

% waist = [623, 1064; 700, 1064]*1e-9;
waist = [623, 1064]*1e-9;

t_move_ramp = [0e-3 0.02e-3];

x = linspace(-0.5e-6,0.5e-6,Nx);
z = 1e-6; %linspace(0,3e-6,Nz);

phi = [7*pi/180 0];
% power = [0.4*1.8e-3 6.8e-3];
power = [10*1.8e-3 6.8e-3];

blast_boo = 1;

[X,Z] = ndgrid(x,z);
xx = X(:);
zz = Z(:);

for j = 1:size(waist,1)
    clear gs_overlap E_out params P_lost;
    
    start_time = now;
    
    for i = 1:numel(zz)
        disp([num2str(i) '/' num2str(numel(zz))])
        separation = [xx(i) zz(i)];
        [gs_overlap(i,:),E_out(i,:),params(i),P_lost(i,:)] = merge_fg_splitstep_2d_fun(separation,waist(j,:),t_move_ramp,power,phi,blast_boo);
        disp(P_lost(i,:))
    end
    params = struct2table(params);
    
    save(['data_' datestr(start_time,'YYmmDD_HHMMSS') '.mat'],'params','gs_overlap','E_out','P_lost')
    
    gs_overlap = reshape(permute(gs_overlap,[1 3 2]),[Nx,Nz,2]);
    
    figure(1);
    clf;
    plot(params.align_err(:,1)*1e6,1-P_lost(:,2),'.-')
    xlabel('radial misalignment (um)')
    ylabel('Cs survival')
    title({['\lambda_{blast} = ' num2str(waist(1)*1e9) ' nm, \phi = ' num2str(phi(1)*180/pi) ' deg'],...
        ['\Deltaz = ' num2str(z*1e6) ' um, P_{blast} = ' num2str(power(1)*1e3) ' mW ']})
    
%     figure(1);
%     clf;
%     subplot(1,2,1);
%     imagesc(z*1e6,x*1e6,gs_overlap(:,:,1));
%     shading flat
%     xlabel('z (um)');
%     xlabel('x (um)');
%     colorbar
%     
%     subplot(1,2,2);
%     imagesc(z*1e6,x*1e6,gs_overlap(:,:,2));
%     shading flat
%     xlabel('z (um)');
%     xlabel('x (um)');
%     colorbar
end