


%% rabi frequencies per sqrt(mW)
switch basis
    case 'aFC'
        [q_ind,v_ind] = evec_ind({'J','I','F','m_F'},[1,5,6,5],c,c.psi);
    case 'aUC'
        [q_ind,v_ind] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[1,1,3/2,5/2],c,c.psi);
    case 'aIC'
        [q_ind,v_ind] = evec_ind({'J','m_J','I','m_I'},[1,1,5,4],c,c.psi);
end
rabi_up_sqrt_mW = max(abs(2*H_up(v_ind)))/const.h*1e-6 *sqrt(1e-3/power_up);
fprintf('up leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_up_sqrt_mW)

switch basis
    case 'aFC'
        [~,v_ind2] = evec_ind({'J','I','F','m_F'},[0,5,5,4],X,X.psi);
    case 'aIC'
        [~,v_ind2] = evec_ind({'J','m_J','I','m_I'},[0,0,5,4],X,X.psi);
    case 'aUC'
        [~,v_ind2] = evec_ind({'J','m_J','m_i_Na','m_i_Cs'},[0,0,3/2,5/2],X,X.psi);
end
rabi_dn_sqrt_mW = max(abs(2*H_dn(v_ind,v_ind2)))/const.h*1e-6 *sqrt(1e-3/power_dn);
fprintf('down leg transition strength = %1.3g MHz/sqrt(mW)\n',rabi_dn_sqrt_mW)


%% list states by raman coupling

Rsc_line =@(Delta,Omega,Gamma) Gamma/2 .* 2*Omega.^2./(2*Omega.^2 + 4*Delta.^2 + Gamma.^2);

raman_amp = abs(H_up.*H_dn);
raman_max = 1.1150e-51; % hard-coded from sig+/sig+ J,M=1,1 to J,M=0,0
raman_amp = raman_amp/raman_max;

[val,order] = sort(raman_amp(:),'descend');

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



%% cut states that aren't coupled
% cut X states that aren't accessed via raman
X_cut = sqrt(abs(H_dn'*H_up))/const.h < 0.01*max(sqrt(abs(H_dn'*H_up))/const.h);
X.E(X_cut) = [];
X.psi(:,X_cut) = [];
H_dn(:,X_cut) = [];

% cut qnums in the ground state that don't have significant admixtures
X_qnums_cut = all(abs(X.psi).^2<1e-3,2);
X.qnums(X_qnums_cut,:) = [];
X.psi(X_qnums_cut,:) = [];
X.psi = X.psi./sqrt(sum(abs(X.psi).^2,1));

% cut c states that don't couple strongly to either ground state
c_cut = (all(abs(H_dn) < 0.1*max(abs(H_dn(:))),2)) & (abs(H_up) < 0.1*max(abs(H_up)));
c.psi(:,c_cut) = [];
c.E(c_cut) = [];
H_up(c_cut,:) = [];
H_dn(c_cut,:) = [];

% cut qnums in the excited state that don't have significant admixtures
c_qnums_cut = all(abs(c.psi).^2<1e-3,2);
c.qnums(c_qnums_cut,:) = [];
c.psi(c_qnums_cut,:) = [];
c.psi = c.psi./sqrt(sum(abs(c.psi).^2,1));

% f.qnums.amp = trapz(f.r,abs(f.psi(:,:,1)).^2,2);

%% build effective hamiltonian
Nf = 1;
Nc = size(c.psi,2);
NX = size(X.psi,2);

Nstates = Nf + Nc + NX;

b1.qnums = build_basis({'eta','i'},{1,1},[0,0]);
b2.qnums = build_basis({'eta','i'},{2,1:size(c.psi,2)},[0,0]);
b3.qnums = build_basis({'eta','i'},{3,1:size(X.psi,2)},[0,0]);
b = join_basis(join_basis(b1,b2),b3);

b.ops.H0 = diag(cat(1,0,c.E-E0_up-min(f.E),X.E-(E0_up-E0_dn)-min(f.E)));
b.ops.Hu = [zeros(Nf) H_up' zeros(Nf,NX); H_up zeros(Nc,Nc+NX); zeros(NX,Nstates)];
b.ops.Hd = [zeros(Nf,Nstates); zeros(Nc,Nf+Nc) H_dn; zeros(NX,Nf) H_dn' zeros(NX)];

b.ops.If = diag(cat(1,ones(Nf,1),zeros(Nc,1),zeros(NX,1)));
b.ops.Ic = diag(cat(1,zeros(Nf,1),ones(Nc,1),zeros(NX,1)));
b.ops.IX = diag(cat(1,zeros(Nf,1),zeros(Nc,1),ones(NX,1)));


%% calculate raman parameters

% detuning_1ph = linspace(-25e9,25e9,1e3);
% 
Gamma = const.c3Sigma.Gamma;
Rsc_bg = 1/100e-6;
% 
% P = b.ops.If + b.ops.IX;
% Q = eye(size(P))-P;
% 
% V = b.ops.Hu+b.ops.Hd;
% 
% ACstark_tot = zeros(size(detuning_1ph));
% RamanRabi = zeros(NX,numel(detuning_1ph));
% Rsc = zeros(size(detuning_1ph));
% tau = zeros(size(detuning_1ph));
% 
% for i = 1:numel(detuning_1ph)
%     H0 = b.ops.H0 - (detuning_1ph(i)*const.h-1i*Gamma)*b.ops.Ic;
%     
%     detuning = H0(b.ops.Ic>0);
%     
%     Gp = zeros(size(H0));
%     Gp(b.ops.Ic>0) = -1./detuning;
%     
% %     Hmm = P*H0*P;
% %     Hpp = Q*H0*Q;
%     
%     Vmm = P*V*P;
%     Vmp = P*V*Q;
%     Vpm = Q*V*P;
%     Vpp = Q*V*Q;
%     
%     Heff = Vmp*Gp*Vpm;
%     
%     ACstark_tot(i) = real(Heff(end,end)-Heff(1,1))/const.h;
%     RamanRabi(:,i) = 2*real(Heff(1,find(diag(b.ops.IX))))/const.hbar;
%     Rsc(i) = (imag(Heff(end,end)+Heff(1,1))/const.hbar) + Rsc_bg;
%     tau(i) = 1./Rsc(i);
% end
% 
% figure(1);
% clf;
% hold on;
% box on;
% plot(detuning_1ph*1e-9,RamanRabi./Rsc,'linewidth',2)
% % plot(detuning_1ph*1e-9,Rsc*1e-8,'-k','linewidth',1)
% hold off;
% grid on
% set(gca,'fontsize',14)
% xlabel('Detuning \Delta (GHz)')
% ylabel('Coherence ratio \Omega_R/R_{sc}')
% xlim([-25 25])
% 
% figure(2);
% clf;
% hold on;
% box on;
% % plot(detuning_1ph*1e-9,RamanRabi./Rsc,'linewidth',2)
% plot(detuning_1ph*1e-9,Rsc*1e-6,'-k','linewidth',1)
% hold off;
% grid on
% set(gca,'fontsize',14)
% xlabel('Detuning \Delta (GHz)')
% ylabel('R_{sc} \times 10^{-6}')
% xlim([-25 25])
% 
% X_qnums_numerator = X.qnums{:,{'m_i_Na','m_i_Cs'}}*2;
% for i = 1:size(X_qnums_numerator,1)
%     temp = num2cell(X_qnums_numerator(i,:));
%     leg{i} = sprintf('m_i_Na = %d/2, m_i_Cs = %d/2',temp{:});
% end
% set(legend(leg,'location','best'),'interpreter','none');
% 
% % fprintf('tau = %0.1f us\nAC stark = %0.1f kHz\nRaman Rabi = %0.1f kHz\n',tau,ACstark_tot,RamanRabi)
% 
% %% calculate raman parameters
% 
% detuning_1ph = 0; % linspace(-25e9,25e9,1e3);
% detuning_2ph = linspace(-25e9,25e9,1e3);
% 
% Gamma = const.c3Sigma.Gamma;
% Rsc_bg = 1/100e-6;
% 
% P = b.ops.If + b.ops.IX;
% Q = eye(size(P))-P;
% 
% V = b.ops.Hu+b.ops.Hd;
% 
% ACstark_tot = zeros(size(detuning_2ph));
% RamanRabi = zeros(NX,numel(detuning_2ph));
% Rsc = zeros(size(detuning_2ph));
% tau = zeros(size(detuning_2ph));
% 
% for i = 1:numel(detuning_2ph)
%     H0 = b.ops.H0 - (detuning_1ph*const.h-1i*Gamma)*b.ops.Ic + const.h*detuning_2ph(i)*b.ops.IX;
%     
%     detuning = H0(b.ops.Ic>0);
%     
%     Gp = zeros(size(H0));
%     Gp(b.ops.Ic>0) = -1./detuning;
%     
% %     Hmm = P*H0*P;
% %     Hpp = Q*H0*Q;
%     
%     Vmm = P*V*P;
%     Vmp = P*V*Q;
%     Vpm = Q*V*P;
%     Vpp = Q*V*Q;
%     
%     Heff = Vmp*Gp*Vpm;
%     
%     ACstark_tot(i) = real(Heff(end,end)-Heff(1,1))/const.h;
%     RamanRabi(:,i) = 2*real(Heff(1,find(diag(b.ops.IX))))/const.hbar;
%     Rsc(i) = (imag(Heff(end,end)+Heff(1,1))/const.hbar) + Rsc_bg;
%     tau(i) = 1./Rsc(i);
% end
% 
% figure(3);
% clf;
% hold on;
% box on;
% plot(detuning_2ph*1e-9,RamanRabi./Rsc,'linewidth',2)
% % plot(detuning_1ph*1e-9,Rsc*1e-8,'-k','linewidth',1)
% hold off;
% grid on
% set(gca,'fontsize',14)
% xlabel('Detuning \Delta (GHz)')
% ylabel('Coherence ratio \Omega_R/R_{sc}')
% xlim([-25 25])
% 
% X_qnums_numerator = X.qnums{:,{'m_i_Na','m_i_Cs'}}*2;
% for i = 1:size(X_qnums_numerator,1)
%     temp = num2cell(X_qnums_numerator(i,:));
%     leg{i} = sprintf('m_i_Na = %d/2, m_i_Cs = %d/2',temp{:});
% end
% set(legend(leg,'location','best'),'interpreter','none');
% 
% figure(4);
% clf;
% hold on;
% box on;
% % plot(detuning_1ph*1e-9,RamanRabi./Rsc,'linewidth',2)
% plot(detuning_2ph*1e-9,Rsc*1e-6,'-k','linewidth',1)
% hold off;
% grid on
% set(gca,'fontsize',14)
% xlabel('Detuning \Delta (GHz)')
% ylabel('R_{sc} \times 10^{-6}')
% xlim([-25 25])
% 
% 
% % fprintf('tau = %0.1f us\nAC stark = %0.1f kHz\nRaman Rabi = %0.1f kHz\n',tau,ACstark_tot,RamanRabi)
