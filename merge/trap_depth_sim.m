clear;

x = linspace(0,0.1e-6,2);

[x,z] = ndgrid(x,z);
x = x(:);
z = z(:);

start_time = now;

for i = 1:numel(z)
    separation = [x(i) z(i)];
    [gs_overlap(i,:),E_out(i,:)] = merge_fg_splitstep_2d_fun(separation);
end

save(['data_' datestr(start_time,'YYMMDD_HHmmSS') '.mat'],'x','z','gs_overlap','E_out')

figure(1);
clf;
% subplot(2,1,1);
plot(sep_z*1e6,gs_overlap,'o-');
legend('Na','Cs');
ylabel('ground state overlap')

% subplot(2,1,2);
% plot(sep_z*1e6,E_out,'o-');
% legend('Na','Cs');
xlabel('z misalignment (um)')
% ylabel('energy (E_z)');
