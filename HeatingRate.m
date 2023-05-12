%everything in SI!
%load constants
C = constants();

auToCmV = 1.648773e-41;
lambda = 1064e-9;
k = 2*pi/lambda;
w0 = lambda*1.17; %
% ItoU_mol = ItoU_Cs*molPol/csPol;

ns = 0:1:20;

T0 = 0; %Initial temp in Kelvin
E0 = 3/2*C.kb*T0;

% I = 2.5e10; %W / m^2. cooling depth
% I = 0.44*2.5e10; %W / m^2. cooling depth

csPol = 1163.4*auToCmV;
% molPol = 939.8*auToCmV;%au

mNaCs = C.m_nacs;
mCs = C.m_cs;

UCsCool = fRadToU(mCs,240e3/2,w0);
UCs = fRadToU(mCs,194.2e3/2,w0);
Umol = fRadToU(mNaCs,171.7e3/2,w0);

I = ItoU(UCs,csPol,1);
IScat = I*1.834/1.5;
IHold = I/15;
ICool = ItoU(UCsCool,csPol,1);
testU = ItoU(I,csPol);

molPol = csPol*(Umol/UCs);

% % MHz/(W/m^2)
%  
% Umol = ItoU_mol*I;
% UCs = ItoU_Cs*I;

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
    xMean(end+1) = x2(20,radFreq,Ts(i),radFreq);
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
    res = sqrt(2*fRad*lambda)/(2*fAx*pi);
end

function res = fRadToU(m,fRad,w0) %Convert radial trap frequency to a trap depth in Hz
    C = constants();
    res = m.*(2*pi.*fRad).^2.*w0.^2/4/C.h; 
end

function res = x2(nMax,trapFreq,T,m)
    C = constants();
    alpha = exp(-C.hbar*2*pi*trapFreq./(C.kb*T));
    pn = sqrt(1 - alpha).*alpha.^(0:nMax);
    xs = sqrt(C.hbar/(m*2*pi*trapFreq).*((0:nMax) + 1/2));
    res = sum(pn.*xs);
end