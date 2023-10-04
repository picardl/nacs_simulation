using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using LsqFit
using Profile
using Random
using DifferentialEquations

N = 2
tPi = 20.5e-6; #Microwave pi pulse time at zero detuning
inter_t = 1.32e-3;

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12
SR = 1e6;

NGauss = 100000;

dets = [500,600]
Δ = 2*pi*500;

gaussStd = 2*pi*10000;
gaussNoise = gaussStd*randn(NGauss)
#gaussFunc = LinearInterpolation(range(0,1,NGauss), gaussNoise)
gaussFunc = (t) -> gaussStd*sin.(2*pi*1.055e3.*t) .+ 1e3
 
#figure(1)
#plot(range(0,50e-3,5000),gaussFunc(range(0,50e-3,5000)))

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

Δt1 = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    #2*pi*[dets[1]]*fracIntens(t); 
    #[gaussFunc(t)];
    [2*pi*500]*(10*fracIntens(t) - 9);
end
Δt2 = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    2*pi*[dets[2]]*fracIntens(t); 
end
Ωt = (t) -> begin
    #2*pi*(1/(4*tPi))*[1,2.805,0.58];
    #2*pi*(1/(4*tPi))*[1,0.2725,0.6585];
    2*pi*(1/(4*tPi))*[1];
end


XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperatorsTimeDep(N, Ωt, Δt1)
#=b = NLevelBasis(N)
XRot1 = 2*pi*(1/(4*tPi)) * (transition(b, 1, 2) + dagger(transition(b, 1, 2)))
YRot1 = 2*pi*(1/(4*tPi)) * (-im*transition(b, 1, 2) + im*dagger(transition(b, 1, 2)))
Ps1 = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:2]
FreeEv1 = Δ * Ps1[2]=#

if true
#Bare Ramsey
psi0 = nlevelstate(b1,1)
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = range(0,2*pi,15);

tWaits = range(0,50e-3,1001)

pf_waits = Array{Any, 1}(undef, length(tWaits));
ψ_wait = Array{Any, 1}(undef, length(tWaits));
pf_stat = Array{Any, 1}(undef, length(tWaits));
ψ_stat = Array{Any, 1}(undef, length(tWaits));

ψfN = Array{Any, 1}(undef, length(tWaits));
pfN = Array{Any, 1}(undef, length(tWaits));
pfNStat = Array{Any, 1}(undef, length(tWaits));

psi0 = nlevelstate(b1,1)
fMod = (t)->(10*fracIntens(t)-9)
HEV = (t,psi) -> FreeEv1*fMod(t);

#=
tout, ψ_t = timeevolution.schroedinger([0,tPi/2], psi0, XRot1 + FreeEv1)
ψ_sp = ψ_t[end]
pf_waits[1:500], ψ_wait[1:500] = timeevolution.schroedinger_dynamic(tWaits[1:500], ψ_sp, HEV)
pf_stat[1:500], ψ_stat[1:500] = timeevolution.schroedinger(tWaits[1:500], ψ_sp, FreeEv1)
tout, ψf = timeevolution.schroedinger([0,tPi], ψ_wait[500], XRot1+ FreeEv1)
toutStat, ψfStat = timeevolution.schroedinger([0,tPi], ψ_stat[500], XRot1+ FreeEv1)
fMod = (t)->(10*fracIntens(t + tWaits[500])-9)
HEV = TimeDependentSum(fMod=>FreeEv1);
pf_waits[500:1001], ψ_wait[500:1001] = timeevolution.schroedinger_dynamic(tWaits[1:502], ψf[end], HEV)
pf_stat[500:1001], ψ_stat[500:1001] = timeevolution.schroedinger(tWaits[1:502], ψfStat[end], FreeEv1)
=#
SEWaits = [23e-3,23e-3]
tout, ψ_wait = DynamicalDecoupling.DDSeqTimeDep([tPi/2,tPi],SEWaits,[0,0],XRot1,YRot1,FreeEv1,psi0,1e6)

pfN = Array{Any, 1}(undef, length(tout));
pfNStat = Array{Any, 1}(undef, length(tout));

H_t = (t,psi) -> begin
    FreeEv1(t,psi) + XRot1(t,psi)
end

for i = 1:length(tout)
    tout2, ψf = timeevolution.schroedinger_dynamic([0,tPi/2], ψ_wait[i],H_t)
    #ψfN[i] = ψf[end];
    pfN[i] = real(expect(Ps1[1],ψf[end]))

    #_, ψf_stat = timeevolution.schroedinger_dynamic([0,tPi/2], ψ_stat[i], XRot1 + FreeEv1)
    #pfNStat[i] = real(expect(Ps1[1],ψf_stat[end]))
end

figure(2)
plot(tout,pfN,label="noisy")
#plot(tWaits,pfNStat,label="fixed")
legend()
xlabel("Ramsey popn")
ylabel("Popn")
ylim([-0.01,1.01])

#=
figure(3)

probePhases = range(0,2*pi,20);
pp = zeros(size(probePhases));
ψp = Array{Any, 1}(undef, length(probePhases))
for i = 1:length(probePhases)
    toutf, psif = timeevolution.schroedinger([0,tPi/2], ψ_wait[end], FreeEv1 + cos(probePhases[i])*XRot1+ sin(probePhases[i])*YRot1)
    pp[i] = real(expect(Ps1[1], psif[end]))
    ψp[i] = psif[end]
end
plot(probePhases,pp)
legend()
xlabel("Ramsey phase")
ylabel("Popn")
ylim([-0.01,1.01])
=#
end

if false

    figure(5)
    tsRabi = range(0,10*tPi,200)
    toutr, psi_t = timeevolution.schroedinger(tsRabi, psi0,XRot1 + FreeEv1)
    exp_val = expect(Ps1[1], psi_t)
    plot(toutr,exp_val,label="|0⟩")
    xlabel("Time (s)")
    ylabel("Popn in |0>")
    title("Rabi")
    legend()

end
