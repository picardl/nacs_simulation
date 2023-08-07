%everything in SI!
%load constants
C = constants();

freqsToWaist(160e3,32e3,1064e-9)

auToCmV = 1.648773e-41; %
lambda = 1064e-9;
k = 2*pi/lambda;
% w0 = lambda*1.17; %Calibrated from axial to radial ratio
w0 = 1.248e-6;%+/-3e-9 From 20230420_105116 param heating

ns = 0:1:20;

T0 = 0; %Initial temp in Kelvin
E0 = 3/2*C.kb*T0;

csPol = 1163.4*auToCmV;

mNaCs = C.m_nacs;
mCs = C.m_cs;

% UCsCool = fRadToU(mCs,218e3/2,1.19e-6); %Cooling depth
UCsCool = fRadToU(mCs,291e3/2,w0); %VODT of 3.4, pow2amp of 5.9982 from 20230509_192841 param heating
UCs = fRadToU(mCs,194.2e3/2,w0); %Calibration depth for molecule polarizability measurement
Umol = fRadToU(mNaCs,171.7e3/2,w0);

UCsAx = fAxToU(mCs,28.2e3,lambda,1.1974e-06);

I = ItoU(UCs,csPol,1);
IScat = I*1.834/1.5; %Intensity at 40 MHz, roughly 1.834 V ODT for scattering
IHold = I/15;
ICool = ItoU(UCsCool,csPol,1);

molPol = csPol*(Umol/UCs);
disp(['NaCs polarizability ',num2str(molPol/auToCmV),' au'])

radFreq = sqrt(4*C.h.*Umol/(mNaCs*w0^2))/(2*pi);
axFreq = sqrt(2*C.h.*Umol*lambda^2/(mNaCs*pi^2*w0^4))/(2*pi);

csFreq = sqrt(4*C.h.*UCs/(mCs*w0^2))/(2*pi);

%Get recoil energy
wr = C.hbar.^2*k^2/(2*mNaCs);

lossRateI = 2.3e-10; %s^-1 / (W m^-2). From our ground state paper
% lossRateU = lossRateI/ItoU_mol; %s^-1 / (Hz)

scatRateI = 54e-10; %s^-1 / (W m^-2). From arxiv.org/pdf/1707.02168.pdf
% scatRateU = scatRateI/ItoU_mol; %s^-1 / (Hz)

ts = [0:1:200]*1e-3;

Es = E0 + wr*IScat*scatRateI*ts;
Es_MHz = Es/1e6/(C.h);
Ts = (2/3)*Es/C.kb;
% plot(ts,Es_MHz);

figure(1)
subplot(1,3,1)
plot(ts,Ts*1e6);
xlabel('time [s]')
ylabel('Temperature [uK]')

nbars = exp(-C.hbar*2*pi*22e3./(C.kb*Ts))./(1-exp(-C.hbar*2*pi*22e3./(C.kb*Ts)));
subplot(1,3,2)
plot(ts,nbars)
xlabel('time [s]')
ylabel('nbar')

xMean = [];
for i = 1:length(ts)
    xMean(end+1) = x2(20,radFreq,Ts(i),mNaCs);
end

figure(1)
subplot(1,3,3)
plot(ts,xMean*1e9)
xlabel('time [s]')
ylabel('sqrt(<x^2>) [nm]')

xMeanLow = [];
for i = 1:length(ts)
    xMeanLow(end+1) = x2(20,radFreq/sqrt(40),Ts(i),mNaCs);
end
figure(4)
plot(ts,xMeanLow*1e9)
xlabel('time [s]')
ylabel('sqrt(<x^2>) [nm]')

function res = ItoU(x,pol,bInv) %Return trap depth in Hz for an intensity in W / m^2
    if nargin < 3
        bInv = 0;
    end
    C = constants();
    if ~bInv
        res = 1/(2*C.eps0*C.c)*pol/C.h*x;
    else
        res = 1/(1/(2*C.eps0*C.c)*pol/(C.h*x));
    end
end

function res = freqsToWaist(fRad,fAx,lambda)
    res = lambda*(fRad/fAx)/(pi*sqrt(2));
end

function res = fRadToU(m,fRad,w0) %Convert radial trap frequency to a trap depth in Hz
    C = constants();
    res = m.*(2*pi.*fRad).^2.*w0.^2/4/C.h; 
end

function res = fAxToU(m,fAx,lambda,w0) %Convert radial trap frequency to a trap depth in Hz
    C = constants();
    res = (2*pi*fAx)^2*m*pi^2*w0^4/(2*lambda^2)/C.h; 
end

function res = x2(nMax,trapFreq,T,m)
    C = constants();
    alpha = exp(-C.hbar*2*pi*trapFreq./(C.kb*T));
    pn = sqrt(1 - alpha).*alpha.^(0:nMax);
    xs = sqrt(C.hbar/(m*2*pi*trapFreq).*((0:nMax) + 1/2));
    res = sum(pn.*xs);
end