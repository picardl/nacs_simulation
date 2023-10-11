using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

N = 4
Δ = 2*pi*[0,-329.6e3,320e3]; #Microwave detuning
tPi = 10e-6;#34e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))*[1,2.805,0.58]

XRot, YRot, FreeEv, Ps, b = genNLevelOperators(N, Ω, Δ)

tspan = range(0,1e-3,10000)

#Diagonal terms
H = FreeEv + XRot #Microwave drive Hamiltonian

psi = nlevelstate(b,1)

#Evolve single particle Hamiltonian Rabi
tout, psi_t = timeevolution.schroedinger(tspan, psi, H)
exp_val_N0 = expect(Ps[1], psi_t)
exp_val_N10 = expect(Ps[2], psi_t)
exp_val_N1m1 = expect(Ps[3], psi_t)
exp_val_N1p1 = expect(Ps[4], psi_t)

figure(1)
plot(tout,exp_val_N0,label="N=0")
plot(tout,exp_val_N10,label="N=1,0")
plot(tout,exp_val_N1m1,label="N=1,-1")
plot(tout,exp_val_N1p1,label="N=1,1")
legend()
xlabel("Time (s)")
ylabel("Popn")
title("Rabi")

#XY8 phase scan
tsXY,tWaitsXY,phasesXY = genXY8(tPi,200e-6,1)
probePhases = range(0,2*pi,30)

pf, ψf = RamseyPhase(probePhases,tPi,tsXY,tWaitsXY,phasesXY,XRot,YRot,FreeEv,Ps[1],psi,1e6)

figure(2)
plot(probePhases,expect(Ps[1], ψf),label="N=0")
plot(probePhases,expect(Ps[2], ψf),label="N=1,0")
plot(probePhases,expect(Ps[3], ψf),label="N=1,-1")
plot(probePhases,expect(Ps[4], ψf),label="N=1,1")

legend()
xlabel("Probe phase")
ylabel("Popn")
ylim([0,1])
title("XY8")

#XY8 N groups scan
NMax = 20
tau = 0.1e-3;
pfN = Array{Any, 1}(undef, NMax);
ψfN = Array{Any, 1}(undef, NMax);
probePhases = pi
for i = 1:1:NMax
    tsXY,tWaitsXY,phasesXY = genXY8(tPi,tau,i)
    pfN[i], ψfN[i] = RamseyPhase(probePhases,tPi,tsXY,tWaitsXY,phasesXY,XRot,YRot,FreeEv,Ps[1],psi,10e6)
end
figure(3)
plot(1:NMax,expect(Ps[1], ψfN),label="N=0")
plot(1:NMax,expect(Ps[2], ψfN),label="N=1,0")
plot(1:NMax,expect(Ps[3], ψfN),label="N=1,-1")
plot(1:NMax,expect(Ps[4], ψfN),label="N=1,1")

legend()
xlabel("N XY8 groups")
ylabel("Popn")
ylim([0,1])
title("XY8")