function stirap_solve

const = constants();

data = load('data/transfer_Heff_aFC_sigp_sigp.mat');
b = data.b;

fp = 1;
power_up = 10e-3*fp; % W
power_dn = 62e-6*fp; % W
waist = 13e-6; % m
detuning_1ph = 0.7e9; % Hz

%%
% electric fields
Efield_up = sqrt(4*const.eta0*power_up/(pi*waist^2)); % V/m
Efield_dn = sqrt(4*const.eta0*power_dn/(pi*waist^2)); % V/m

%% further eliminate off-resonant states
H0 = b.ops.H0 - const.h*detuning_1ph*b.ops.Ic - 1i*const.c3Sigma.Gamma*b.ops.Ic;
Hu_unit = b.ops.Hu/2;
Hd_unit = b.ops.Hd/2;

Nt = 2e3;
f = 2;
tmax = 100e-6 / f;
thold = 40e-6 / f;
td = 8e-6 / f;
pw = 10e-6 / f;
detuning_2ph = 0; %-0.7e6;

pulse =@(t0,pw,t) exp(-(t-t0).^2./(2*(pw/(2*sqrt(2*log(2))))^2));

E_up =@(t) Efield_up*(pulse(tmax/2-thold/2+td/2,pw,t) + pulse(tmax/2+thold/2-td/2,pw,t));
E_dn =@(t) Efield_dn*(pulse(tmax/2-thold/2-td/2,pw,t) + pulse(tmax/2+thold/2+td/2,pw,t));

t = reshape(linspace(0,tmax,Nt),1,1,[]);
dt = t(2)-t(1);

H = (H0 + Hu_unit.*E_up(t) + Hd_unit.*E_dn(t) - const.h*detuning_2ph*b.ops.IX)/const.hbar;
Nstates = size(H,1);

prop = cell2mat(arrayfun(@(i) expm(-1i*H(:,:,i)*dt),reshape(1:Nt-1,1,1,Nt-1),'un',0));

psi = zeros(Nstates,Nt);
psi(1,1) = 1;
for i = 2:Nt
    psi(:,i) = prop(:,:,i-1)*psi(:,i-1);
end

figure(1);
clf;
subplot(2,1,1);
hold on;
box on;
plot(t(:),E_up(t(:))/Efield_up);
plot(t(:),E_dn(t(:))/Efield_dn);
hold off;
legend('pump','stokes')

subplot(2,1,2);
plot(t(:),abs(psi).^2);
ylim([0 1])
legend(num2str(b.qnums{:,:}))

%%


end