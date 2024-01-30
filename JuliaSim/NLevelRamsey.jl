using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

N = 4
U = 3.5e6;
Δ = 2*pi*[0,-0.2*U,0.2*U]; #Microwave detuning
tPi = 10.202e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))*[1,0.2724,0.658] #Based on dataset 20230811_101049

XRot, YRot, FreeEv, Ps, b = genNLevelOperators(N, Ω, Δ)

tspan = range(0,1e-3,10000)

#Use general pulse sequence formulation for spin echo

if false
    psi = nlevelstate(b,1)

    amps = [1,0,1,0,1];
    times = tPi*[1/2,10,1,10,1/2];
    phases::Vector{Float64} = [0,0,0,0,0];

    params = Dict();
    params["frac"] = 0.3;
    probePhases = range(0,2*pi,30)
    ψf = Array{Any, 1}(undef, length(probePhases));

    for i = 1:length(probePhases)
        phases[end] = probePhases[i];
        tTots,ψTots,Ps = DynamicalDecoupling.generalPulseSeq(psi,amps,times,phases,Ω,Δ,N,"square",params);
        ψf[i] = ψTots[end];
    end


    figure(2)
    plot(probePhases,expect(Ps[1], ψf),label="N=0")
    plot(probePhases,expect(Ps[2], ψf),label="N=1,0")
    plot(probePhases,expect(Ps[3], ψf),label="N=1,-1")
    plot(probePhases,expect(Ps[4], ψf),label="N=1,1")

    legend()
    xlabel("Probe phase")
    ylabel("Popn")
    ylim([0,1])
    title("Ramsey spin echo")
end

#Use general pulse sequence formulation for XY8
if false
psi = nlevelstate(b,1)
tTau = 84e-6;

amps_outer = [1,1];
times_outer = tPi*[1/2,1/2];
phases_outer::Vector{Float64} = [0,0];
params = Dict();
params["frac"] = 0.3;

amps_xy8_core = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
times_xy8_core= [tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,tTau]
phases_xy8_core = pi*[0,0,0,0.5,0,0,0,0.5,0,0.5,0,0,0,0.5,0,0,0];

nGroups = 1:5:20;
ψf = Array{Any, 1}(undef, length(nGroups));

for i = 1:length(nGroups)
    print('.')
    amps = vcat(amps_outer[1],repeat(amps_xy8_core,nGroups[i]),amps_outer[2]);
    times = vcat(times_outer[1],repeat(times_xy8_core,nGroups[i]),times_outer[2]);
    phases = vcat(phases_outer[1],repeat(phases_xy8_core,nGroups[i]),phases_outer[2]);
    tTots,ψTots,Ps = DynamicalDecoupling.generalPulseSeq(psi,amps,times,phases,Ω,Δ,N,"square",params);
    ψf[i] = ψTots[end];
end

figure(3)
    plot(nGroups,expect(Ps[1], ψf),label="N=0")
    plot(nGroups,expect(Ps[2], ψf),label="N=1,0")
    plot(nGroups,expect(Ps[3], ψf),label="N=1,-1")
    plot(nGroups,expect(Ps[4], ψf),label="N=1,1")

    legend()
    xlabel("N Groups")
    ylabel("Popn")
    ylim([0,1])
    title("XY8")
end

#Use stochastic general pulse sequence formulation for XY8
if false
    psi = nlevelstate(b,1)
    tTau = 84e-6;
    
    amps_outer = [1,1];
    times_outer = tPi*[1/2,1/2];
    phases_outer::Vector{Float64} = [0,0];
    
    amps_xy8_core = [0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0];
    times_xy8_core= [tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,tTau]
    phases_xy8_core = pi*[0,0,0,0.5,0,0,0,0.5,0,0.5,0,0,0,0.5,0,0,0];
    
    nGroups = vcat(1:2:20,25:5:80);
    ψf = Array{Any, 1}(undef, length(nGroups));
    params = Dict();
    params["frac"] = 0.3;

    H_noise = FreeEv*0.01;
    
    for i = 1:length(nGroups)
        print('.')
        amps = vcat(amps_outer[1],repeat(amps_xy8_core,nGroups[i]),amps_outer[2]);
        times = vcat(times_outer[1],repeat(times_xy8_core,nGroups[i]),times_outer[2]);
        phases = vcat(phases_outer[1],repeat(phases_xy8_core,nGroups[i]),phases_outer[2]);
        tTots,ψTots,Ps = DynamicalDecoupling.generalPulseSeq_stochastic(psi,amps,times,phases,Ω,Δ,H_noise,N,"truncGaussDiscrete",params);
        ψf[i] = ψTots[end];
    end
    
    figure(4)
        plot(nGroups,expect(Ps[1], ψf),label="N=0")
        plot(nGroups,expect(Ps[2], ψf),label="N=1,0")
        plot(nGroups,expect(Ps[3], ψf),label="N=1,-1")
        plot(nGroups,expect(Ps[4], ψf),label="N=1,1")
    
        legend()
        xlabel("N Groups")
        ylabel("Popn")
        ylim([0,1])
        title("XY8")
end

#Use stochastic general pulse sequence formulation for XY4
if true
    psi = nlevelstate(b,1)
    tTau = 84e-6;
    
    amps_outer = [1,1];
    times_outer = tPi*[1/2,1/2];
    phases_outer::Vector{Float64} = [0,0];
    
    amps_xy4_core = [0,1,0,1,0,1,0,1,0];
    times_xy4_core= [tTau,tPi,2*tTau,tPi,2*tTau,tPi,2*tTau,tPi,tTau]
    phases_xy4_core = pi*[0,0,0,0.5,0,0,0,0.5,0];
    
    nGroups = vcat(1:2:20,25:5:80);
    ψf = Array{Any, 1}(undef, length(nGroups));
    params = Dict();
    params["frac"] = 0.3;

    H_noise = FreeEv*0.01;
    
    for i = 1:length(nGroups)
        print('.')
        amps = vcat(amps_outer[1],repeat(amps_xy4_core,nGroups[i]),amps_outer[2]);
        times = vcat(times_outer[1],repeat(times_xy4_core,nGroups[i]),times_outer[2]);
        phases = vcat(phases_outer[1],repeat(phases_xy4_core,nGroups[i]),phases_outer[2]);
        tTots,ψTots,Ps = DynamicalDecoupling.generalPulseSeq_stochastic(psi,amps,times,phases,Ω,Δ,H_noise,N,"truncGaussDiscrete",params);
        ψf[i] = ψTots[end];
    end
    
    figure(4)
        plot(nGroups,expect(Ps[1], ψf),label="N=0")
        plot(nGroups,expect(Ps[2], ψf),label="N=1,0")
        plot(nGroups,expect(Ps[3], ψf),label="N=1,-1")
        plot(nGroups,expect(Ps[4], ψf),label="N=1,1")
    
        legend()
        xlabel("N Groups")
        ylabel("Popn")
        ylim([0,1])
        title("XY8")
end
#Evolve single particle Hamiltonian Rabi
#=
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
for i = 1:5:NMax
    tsXY,tWaitsXY,phasesXY = genXY8(tPi,200e-6,i)
    pfN[i], ψfN[i] = RamseyPhase(probePhases,tPi,tsXY,tWaitsXY,phasesXY,XRot,YRot,FreeEv,Ps[1],psi,1e6)
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
title("XY8")=#