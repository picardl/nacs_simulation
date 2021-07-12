function raman_spectrum(B)

const = constants;

if nargin<1
    B = 854e-4;
end

raman_data = raman_effective_hamiltonian(B);


%% laser stuff
power_up = 1e-3;
power_dn = 1e-3;

waist = 10e-6;

pol_up = sphten([1 -1i 0]/sqrt(2));
pol_dn = sphten([1 -1i 0]/sqrt(2));

Efield_up = sqrt(4*const.eta0*power_up/(pi*waist^2));
Efield_dn = sqrt(4*const.eta0*power_dn/(pi*waist^2));

%% list states by raman coupling
H_up = Efield_up/2*sum(raman_data.H_up(:,1,:,:).*reshape(pol_up,1,1,1,3),4);
H_dn = Efield_dn/2*sum(raman_data.H_dn.*reshape(pol_dn,1,1,1,1,3),5);

c = raman_data.c;
X = raman_data.X;
f = raman_data.f;

E0_up = min(c.E) - min(f.E);
E0_dn = min(c.E) - min(X.E);

Gamma = const.c3Sigma.Gamma;

Rsc_line =@(Delta,Omega,Gamma) Gamma/2 .* 2*abs(Omega).^2./(2*abs(Omega).^2 + 4*Delta.^2 + Gamma.^2);

raman_amp = abs(H_up.*H_dn);
raman_max = 1.1150e-51; % hard-coded from sig+/sig+ J,M=1,1 to J,M=0,0
raman_amp = raman_amp/raman_max;

[~,order] = sort(raman_amp(:),'descend');

[row,col] = ind2sub(size(raman_amp),order);

raman_amp_sorted = arrayfun(@(i) raman_amp(row(i),col(i)),(1:numel(row))');

up_leg_rabi = abs(H_up(row)/const.h)*1e-6;
down_leg_rabi = abs(arrayfun(@(i) H_dn(row(i),col(i)),(1:numel(row))')/const.h)*1e-6;

upper_qnums = c.qnums(evec_leading_percentages(abs(c.psi(:,row)).^2,1),:);
lower_qnums = X.qnums(evec_leading_percentages(abs(X.psi(:,col)).^2,1),:);

upper_qnums.Properties.VariableNames = strcat(upper_qnums.Properties.VariableNames,'_i');
lower_qnums.Properties.VariableNames = strcat(lower_qnums.Properties.VariableNames,'_f');

detuning_up_rel = 1e-9 * (c.E(row)-min(f.E)-E0_up)/const.h;
detuning_dn_rel = (c.E(row)-X.E(col) - E0_dn)/const.h*1e-9;




detuning_up = (c.E' - c.E(row));

H_up_temp = repmat(H_up(row),1,size(c.psi,2));
H_up_temp(detuning_up==0) = 0;
Rsc_offres_up = 1e-6*sum(Rsc_line(detuning_up,H_up_temp,Gamma),2)/const.h;

H_up_temp = repmat(H_up(row),1,size(c.psi,2));
H_up_temp(detuning_up~=0) = 0;
Rsc_res_up = 1e-6*sum(Rsc_line(detuning_up,H_up_temp,Gamma),2)/const.h;

detuning_dn = (X.E' - X.E(col));

H_dn_temp = arrayfun(@(i) H_dn(row(i),col(i)),(1:numel(row))');
H_dn_temp = repmat(H_dn_temp,1,size(X.psi,2));
H_dn_temp(detuning_dn==0) = 0;
Rsc_offres_dn = 1e-6*sum(Rsc_line(detuning_dn,H_dn_temp,Gamma),2)/const.h;

H_dn_temp = arrayfun(@(i) H_dn(row(i),col(i)),(1:numel(row))');
H_dn_temp = repmat(H_dn_temp,1,size(X.psi,2));
H_dn_temp(detuning_dn~=0) = 0;
Rsc_res_dn = 1e-6*sum(Rsc_line(detuning_dn,H_dn_temp,Gamma),2)/const.h;


results_table = upper_qnums(:,{'J_i','m_J_i','m_i_Na_i','m_i_Cs_i'});
results_table = cat(2,results_table,lower_qnums(:,{'J_f','m_J_f','m_i_Na_f','m_i_Cs_f'}));
results_table.up_rabi = round(up_leg_rabi,4);
results_table.dn_rabi = round(down_leg_rabi,4);
results_table.detuning_up = round(detuning_up_rel,4);
results_table.detuning_dn = round(detuning_dn_rel,4);
results_table.raman_amp = round(raman_amp_sorted,4);
results_table.Rsc_res_up = round(Rsc_res_up,4);
results_table.Rsc_offres_up = round(Rsc_offres_up,4);
results_table.Rsc_res_dn = round(Rsc_res_dn,4);
results_table.Rsc_offres_dn = round(Rsc_offres_dn,4);

disp(results_table(1:10,:))



end
