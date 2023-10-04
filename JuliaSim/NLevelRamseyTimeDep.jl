using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using LsqFit
using Profile

N = 4
tPi = 7.4235e-6; #Microwave pi pulse time at zero detuning

#csv_file_path = "C:/nilab-projects/nacs_simulation/JuliaSim/dcNoise250ms.CSV"
csv_file_path = "C:/projects/nacs_simulation/JuliaSim/dcNoise250ms.CSV"
lines = readlines(csv_file_path)
lines = map(strip, lines)
data = map(line -> split(line, ','), lines)
# Convert the 4th and 5th columns to Float64, leave other columns as strings
data = map(row -> [row[1], row[2], row[3], parse(Float64, row[4]), parse(Float64, row[5]), row[6:end]], data)
# Extract the 4th column
data_t = map(row -> row[4], data)
# Extract the 5th column
data_v = map(row -> row[5], data)
data_t = data_t[100:end] .- data_t[100]
data_v = data_v[100:end]

fracIntens = LinearInterpolation(data_t, data_v./mean(data_v))

Δt = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    2*pi*[0,-261.17e3,261.17e3]*fracIntens(t); 
end
Ωt = (t) -> begin
    #2*pi*(1/(4*tPi))*[1,2.805,0.58];
    #2*pi*(1/(4*tPi))*[1,0.2725,0.6585];
    2*pi*(1/(4*tPi))*[1,1,1];
end

XRot, YRot, FreeEv, Ps, b = genNLevelOperatorsTimeDep(N, Ωt, Δt)

tspan = range(0,1e-3,1000);

#Diagonal terms
Ht = (t,psi) -> begin
    FreeEv(t,psi) + XRot(t,psi) #Microwave drive Hamiltonian
end

psi = nlevelstate(b,1);

#Evolve single particle Hamiltonian Rabi
if true
tout, psi_t = timeevolution.schroedinger_dynamic(tspan, psi, Ht)
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
end
#=
#Spin echo time scan
tWaits = range(0,30e-3,20)
tsSpinEcho = tPi.*tFracSpinEcho
probePhases = 0
pfN = Array{Any, 1}(undef, length(tWaits));
ψfN = Array{Any, 1}(undef, length(tWaits));

for i = 1:1:length(tWaits)
    tWaitsSpinEcho = tWaits[i]*waitFracSpinEcho
    pfN[i], ψfN[i] = RamseyPhaseTimeDep(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,XRot,YRot,FreeEv,Ps[1],psi,1e6)
end

figure(3)
plot(tWaits,expect(Ps[1], ψfN),label="N=0")
plot(tWaits,expect(Ps[2], ψfN),label="N=1,0")
plot(tWaits,expect(Ps[3], ψfN),label="N=1,-1")
plot(tWaits,expect(Ps[4], ψfN),label="N=1,1")

legend()
xlabel("Spin echo wait time")
ylabel("Popn")
ylim([0,1])
=#

#XY8 N groups scan
if false
NMax = 20;
dN = 2;
tau = 0.1e-3;
pfN = Array{Any, 1}(undef, length(1:dN:NMax));
ψfN = Array{Any, 1}(undef, length(1:dN:NMax));
probePhases = pi
global i = 1
for n = 1:dN:NMax
    tsXY,tWaitsXY,phasesXY = DynamicalDecoupling.genXY8(tPi,tau,n)
    pfN[i], ψfN[i] = RamseyPhaseTimeDep(probePhases,tPi,tsXY,tWaitsXY,phasesXY,XRot,YRot,FreeEv,Ps[1],psi,1e6)
    global i = i+1;
end
figure(3)
plot(1:dN:NMax,expect(Ps[1], ψfN),label="N=0")
plot(1:dN:NMax,expect(Ps[2], ψfN),label="N=1,0")
plot(1:dN:NMax,expect(Ps[3], ψfN),label="N=1,-1")
plot(1:dN:NMax,expect(Ps[4], ψfN),label="N=1,1")

legend()
xlabel("N XY8 groups")
ylabel("Popn")
ylim([0,1])
title("XY8")
end

#XY8 N groups scan monte carlo
if true
    NTrials = 10;

    NMax = 20;
    dN = 5;
    tau = 1.25e-3;
    pfN = Array{Any, 1}(undef, length(1:dN:NMax));
    mean_values = zeros(Float64,length(1:dN:NMax),N)
    stderr_values = zeros(Float64,length(1:dN:NMax),N)

    these_expect = zeros(Float64,N,NTrials)

    probePhases = pi
    global i = 1

    tStart = rand(1,NTrials)
    Xoffs = Vector{Any}(undef, NTrials)
    Yoffs = Vector{Any}(undef, NTrials)
    Freeoffs = Vector{Any}(undef, NTrials)
    for j=1:NTrials
        t_offset = tStart[j]
        Xoffs[j] = (t, psi) -> XRot(t + t_offset, psi)
        Yoffs[j] = (t,psi)-> YRot(t + t_offset, psi)
        Freeoffs[j] = (t,psi)-> FreeEv(t + t_offset, psi)
    end

    for n = 1:dN:NMax
        tsXY,tWaitsXY,phasesXY = DynamicalDecoupling.genXY8(tPi,tau,n)
        for j=1:NTrials
            @time pf,ψf = RamseyPhaseTimeDep(probePhases,tPi,tsXY,tWaitsXY,phasesXY,Xoffs[j],Yoffs[j],Freeoffs[j],Ps[1],psi,5e6)
            these_expect[:,j] = [inner[1] for inner in real.(expect.(Ps, Ref(ψf)))]
        end
        mean_values[i,:] = mean(these_expect, dims=2)
        stderr_values[i,:] = std(these_expect, dims=2)./sqrt(NTrials)
        global i = i+1;
    end
    figure(3)
    errorbar(1:dN:NMax,mean_values[:,1], stderr_values[:,1],label="N=0")
    errorbar(1:dN:NMax,mean_values[:,2], stderr_values[:,2],label="N=1,0")
    errorbar(1:dN:NMax,mean_values[:,3], stderr_values[:,3],label="N=1,-1")
    errorbar(1:dN:NMax,mean_values[:,4], stderr_values[:,4],label="N=1,1")
    
    legend()
    xlabel("N XY8 groups")
    ylabel("Popn")
    ylim([0,1])
    title("XY8")
end

#XY8 scan pi time
if false
    NTrials = 20;

    NMax = 30;
    tau = 0.3e-3;

    #piTimes = [1,2,5,10,15,20,25,30,35,40,50]*1e-6;
    piTimes = [1,2,10,20,25,30,40,50,100]*1e-6;

    pfN = Array{Any, 1}(undef, length(piTimes));
    mean_values = zeros(Float64,length(piTimes),N)
    stderr_values = zeros(Float64,length(piTimes),N)

    these_expect = zeros(Float64,N,NTrials)

    probePhases = pi
    global i = 1

    tStart = rand(1,NTrials)
    Xoffs = Vector{Any}(undef, NTrials)
    Yoffs = Vector{Any}(undef, NTrials)
    Freeoffs = Vector{Any}(undef, NTrials)
    for j=1:NTrials

    end

    for n = 1:1:length(piTimes)
        tFrac = piTimes[n]/tPi
        tsXY,tWaitsXY,phasesXY = DynamicalDecoupling.genXY8(tPi,tau,NMax)
        @time begin 
            for j=1:NTrials
                t_offset = tStart[j]
                Xoffs = (t, psi) -> XRot(t + t_offset, psi)/tFrac
                Yoffs = (t,psi)-> YRot(t + t_offset, psi)/tFrac
                Freeoffs = (t,psi)-> FreeEv(t + t_offset, psi)
                pf,ψf = RamseyPhaseTimeDep(probePhases,tPi*tFrac,tsXY*tFrac,tWaitsXY,phasesXY,Xoffs,Yoffs,Freeoffs,Ps[1],psi,5e6)
                these_expect[:,j] = [inner[1] for inner in real.(expect.(Ps, Ref(ψf)))]
            end
        end
        mean_values[i,:] = mean(these_expect, dims=2)
        stderr_values[i,:] = std(these_expect, dims=2)./sqrt(NTrials)
        global i = i+1;
    end
    figure(3)
    errorbar(piTimes,mean_values[:,1], stderr_values[:,1],label="N=0")
    errorbar(piTimes,mean_values[:,2], stderr_values[:,2],label="N=1,0")
    errorbar(piTimes,mean_values[:,3], stderr_values[:,3],label="N=1,-1")
    errorbar(piTimes,mean_values[:,4], stderr_values[:,4],label="N=1,1")
    
    legend()
    xlabel("Pi time")
    ylabel("Popn")
    ylim([0,1])
    title("XY8")
end