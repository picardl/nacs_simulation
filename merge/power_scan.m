clear;

waist = [623, 1064; 700, 1064]*1e-9;
% waist = [700, 1064]*1e-9;

t_move_ramp = [0.1e-3 0.1e-3];

pow = linspace(0,1,31);

for j = 1:size(waist,1)
    clear gs_overlap E_out params;
    
    start_time = now;
    
    for i = 1:numel(pow)
        disp([num2str(i) '/' num2str(numel(pow))])
        [gs_overlap(i,:),E_out(i,:),params(i)] = merge_fg_splitstep_2d_fun([0 0],waist(j,:),t_move_ramp,[pow(i)*1.8e-3 6.8e-3]);
    end
    params = struct2table(params);
    
    save(['data_' datestr(start_time,'YYmmDD_HHMMSS') '.mat'],'params','gs_overlap','E_out')
    
end