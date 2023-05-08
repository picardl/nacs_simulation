%Testbed for rf phases

t = linspace(0,10e-3,10e6+1);
fs = 1/(t(2)-t(1));

freqs = [105:2.5:182.5]*1e6;
ns = 42:42+length(freqs);

N = length(freqs);
amps = ones(size(freqs))/N;
phases2 = 2*pi*rand(1,N);
phases1 = zeros(1,N);
phases3 = -pi*ns.^2/N;

sig1 = sumSig(t,freqs,amps,phases1,0.02);
sig2 = sumSig(t,freqs,amps,phases2,0.02);
sig3 = sumSig(t,freqs,amps,phases3,0.02);
% sig = sig + sig.^2;

[psd1,f] = pspectrum(sig1,t,'FrequencyResolution',0.05e6);
[psd2,~] = pspectrum(sig2,t,'FrequencyResolution',0.05e6);
[psd3,~] = pspectrum(sig3,t,'FrequencyResolution',0.05e6);

psd1 = psd1/max(psd1);
psd2 = psd2/max(psd2);
psd3 = psd3/max(psd3);

peaks = find(ismembertol(f,freqs));

%%
figure(1)
clf
plot(f,psd1);
xlim([100e6,200e6])
hold on
plot(f,psd2);
plot(f,psd3);
legend({'zero shift','random phases','Schroeder phases'});

figure(2)
clf
plot(f(peaks),psd1(peaks));
xlim([100e6,200e6])
hold on
plot(f(peaks),psd2(peaks));
plot(f(peaks),psd3(peaks));
legend({'zero shift','random phases','Schroeder phases'});

function out = sumSig(t,freqs,amps,phases,nonlin)
    out = 0;
    for i = 1:length(freqs)
        out = out + amps(i).*sin(2*pi*freqs(i)*t + phases(i));
    end
    for i = 1:length(freqs)
        out = out.*(1 + nonlin*sin(2*pi*freqs(i)*t + phases(i)));
    end
end