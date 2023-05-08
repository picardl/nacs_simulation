function upleg_spectrum(B)

const = constants;

if nargin<1
    B = 800e-4;
end

basis = 'aFC';

recompute = 1;

raman_data = raman_effective_hamiltonian(B,basis,recompute,3,[3 4 5],320010);

%% initial state
psi_init = zeros(size(raman_data.f.E));
psi_init(1) = 1;

%% laser stuff
power = 100e-3;
waist = 50e-6;
pol = sphten([1 -1i 0]/sqrt(2));
t = 10e-3;

Efield = sqrt(4*const.eta0*power/(pi*waist^2));

Ndet = 1e3;
freq_offs = 320010e9;
detuning = linspace(0,20e9,Ndet);

%%
rabi = Efield*sum(raman_data.H_up.*psi_init.*permute(pol,[1 3 4 2]),4);

delta = raman_data.c.E - raman_data.f.E - const.h*reshape(freq_offs + detuning,1,1,[]);

Rsc = permute(sum(sum((const.c3Sigma.Gamma/(2*const.h)) * 2*abs(rabi).^2./(2*abs(rabi).^2 + 4*delta.^2 + const.c3Sigma.Gamma.^2),2),1),[2 3 1]);

survival = exp(-Rsc*t);

figure(1); clf;
subplot(2,1,1);
plot(detuning*1e-9,Rsc*1e-6)
xlabel('up leg freq (GHz)')
ylabel('R_{sc} (10^6 s^-1)')

subplot(2,1,2);
plot(detuning*1e-9,survival)
ylim([0 1])
xlabel('up leg freq (GHz)')
ylabel('survival')

end
