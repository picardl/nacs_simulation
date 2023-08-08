using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

#=
function DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    #Does time evolution for a generic dynamical decoupling pulse sequence
    #ts: (N,1) array of rotation times 
    #tWaits: (N,1) array of free evolution times after each pulse in ts
    #phases: (N,1) array of phases of each rotation
    #XRot: Operator representing rotation about X axis (rabi freq should already by factor in)
    #YRot: Operator representing rotation about Y axis (rabi freq should already by factor in)
    #freeEv: Operator representing free evolution between pulses, which should also capture detuning of states during pulse
    #ψ0: initial state Vector
    #SR: sample rate for time points

    tTots = []
    ψTots = []
    ψ = ψ0
    tEnd = 0

    T = 1/SR

    for i = 1:length(ts)
        tPulse = 0:T:ts[i]
        tout, ψ_t = timeevolution.schroedinger(tPulse, ψ, freeEv + cos(phases[i])*XRot + sin(phases[i])*YRot)
        tTots = vcat(tTots,tout .+ tEnd)
        ψTots = vcat(ψTots,ψ_t)

        tEnd = tTots[end]
        ψ = ψTots[end]

        if tWaits[i] > 0
            tPulse = 0:T:tWaits[i]
            tout, ψ_t = timeevolution.schroedinger(tPulse, ψ, freeEv)
            tTots = vcat(tTots,tout .+ tEnd)
            ψTots = vcat(ψTots,ψ_t)

            tEnd = tTots[end]
            ψ = ψTots[end]
        end
    end

    return tTots, ψTots
end

function RamseyPhase(probePhases,tPi,ts,tWaits,phases,XRot,YRot,freeEv,P1,ψ0,SR)
    #Return final |0⟩ popn for final pi/2 phases in probePhases
    tout, psi_t = DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    exp_val = expect(P1, psi_t)
    
    pf = zeros(size(probePhases));
    for i = 1:length(probePhases)
        toutf, psif = timeevolution.schroedinger([0,tPi/2], psi_t[end], freeEv + cos(probePhases[i])*XRot + sin(probePhases[i])*YRot)
        pf[i] = real(expect(P1, psif[end]))
    end
    return pf
end
=#

N = 2
b = NLevelBasis(N)
Δ = 2*pi*(500); #Microwave detuning
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

psi = nlevelstate(b,1)

tspan = range(0,100e-6,10000)

tpi2 =range(0,tPi/2,1000)

#Sigma x and y operators for the two levels
#=
σx = (transition(b,1,2) + dagger(transition(b,1,2)))
σy = (-im*transition(b,1,2) + im*dagger(transition(b,1,2)))

#Projection operators for each state
P1 = tensor(nlevelstate(b,1), dagger(nlevelstate(b,1)))
P2 =  tensor(nlevelstate(b,2), dagger(nlevelstate(b,2)))

#general rotation Hamiltonian
genRot(Ωarg,ϕ) = Ωarg*(cos(ϕ)*σx + sin(ϕ)*σy)
=#
# SIMULATE SPIN ECHO FOR A SINGLE MOLECULE

XRot, YRot, FreeEv, Ps = genNLevelOperators(2, Ω, Δ)
P1 = Ps[1]
probePhases = range(0,2*pi,30)

tsSpinEcho = tPi.*tFracSpinEcho
tWaitsSpinEcho = 360*1e-6*waitFracSpinEcho

pf, ψf = RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,XRot,YRot,FreeEv,P1,psi,1e6)
figure(1)
plot(probePhases,pf)
xlabel("Ramsey phase")
ylabel("N=0 popn")
title("Single molecule spin Echo Ramsey") 


# SIMULATE AVG SPIN ECHO FOR 8 MOLECULES
#=
rabi_frequencies = range(0.7,1.3,8) #Fractional errors in pi time for each state
dets = 500*(rand(8).-0.5).*2*pi;

Xs(rabRat) = Ω*σx*rabRat
Ys(rabRat) = Ω*σy*rabRat
FreeEv(Δ) = Δ*P2

X_col = Xs.(rabi_frequencies)
Y_col = Ys.(rabi_frequencies)
free_col = frees.(dets)

pf_col, allRes = DynamicalDecoupling.CollectiveRamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,X_col,Y_col,free_col,P1,psi,1e6)
figure(2)
plot(probePhases,pf_col)
xlabel("Ramsey phase")
ylabel("N=0 popn")
title("Avg Spinecho Ramsey with 1% Ω err and +- 500 Hz Δ err")=#

# SIMULATE COLLECTIVE SPIN ECHO FOR 8 MOLECULES
#= num_molecules = 8;
rabi_frequencies = range(0.995,1.005,8) #Fractional errors in pi time for each state
dets = 500*(rand(8).-0.5).*2*pi;

b_coll = tensor([b for i=1:num_molecules]...)
Xs(i,rabRat) = embed(b_coll,i,Ω*σx*rabRat)
Ys(i,rabRat) = embed(b_coll,i,Ω*σy*rabRat)
frees(i,Δ) = embed(b_coll,i,Δ*P2)
P1s(i) = embed(b_coll,i,P1)

X_col = sum(Xs.(1:num_molecules,rabi_frequencies))
Y_col = sum(Ys.(1:num_molecules,rabi_frequencies))
free_col = sum(frees.(1:num_molecules,dets))
P_col = sum(P1s.(1:num_molecules))/num_molecules

psi_col = tensor([psi for i=1:num_molecules]...)

pf_col = RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,X_col,Y_col,free_col,P_col,psi_col,1e6)
figure(2)
plot(probePhases,pf_col)
xlabel("Ramsey phase")
ylabel("N=0 popn")
title("Avg Spinecho Ramsey with 1% Ω err and +- 500 Hz Δ err") =#

#= # SIMULATE COLLECTIVE SPIN ECHO FOR 8 MOLECULES DETUNING ONLY
num_molecules = 8;
rabi_frequencies =ones(8) #Fractional errors in pi time for each state
Δ = 0:50:5000

b_coll = tensor([b for i=1:num_molecules]...)
Xs(i,rabRat) = embed(b_coll,i,Ω*σx*rabRat)
Ys(i,rabRat) = embed(b_coll,i,Ω*σy*rabRat)
frees(i,Δ) = embed(b_coll,i,Δ*P2)
P1s(i) = embed(b_coll,i,P1)
P_col = sum(P1s.(1:num_molecules))/num_molecules

X_col = sum(Xs.(1:num_molecules,rabi_frequencies))
Y_col = sum(Ys.(1:num_molecules,rabi_frequencies))

pf_col = zeros(length(Δ))

for i = 1:length(Δ)
    dets = Δ[i]*(range(-0.5,0.5,num_molecules)).*2*pi;
    free_col = sum(frees.(1:num_molecules,dets))

    psi_col = tensor([psi for i=1:num_molecules]...)
    pf_col_loop = RamseyPhase(0,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,X_col,Y_col,free_col,P_col,psi_col,1e6)
    pf_col[i] = pf_col_loop[1]
end
figure(2)
plot(Δ,pf_col)
xlabel("Pk2Pk detuning")
ylabel("Spin echo at π phase") =#

# SIMULATE COLLECTIVE RAMSEY NO ECHO FOR 8 MOLECULES DETUNING ONLY
#=
num_molecules = 8;
rabi_frequencies =ones(8) #Fractional errors in pi time for each state
Δ = 0:10:1500

b_coll = tensor([b for i=1:num_molecules]...)
Xs(i,rabRat) = embed(b_coll,i,Ω*σx*rabRat)
Ys(i,rabRat) = embed(b_coll,i,Ω*σy*rabRat)
frees(i,Δ) = embed(b_coll,i,Δ*P2)
P1s(i) = embed(b_coll,i,P1)
P_col = sum(P1s.(1:num_molecules))/num_molecules

X_col = sum(Xs.(1:num_molecules,rabi_frequencies))
Y_col = sum(Ys.(1:num_molecules,rabi_frequencies))

pf_col = zeros(length(Δ))

for i = 1:length(Δ)
    dets = Δ[i]*(range(-0.5,0.5,num_molecules)).*2*pi;
    free_col = sum(frees.(1:num_molecules,dets))

    psi_col = tensor([psi for i=1:num_molecules]...)
    pf_col_loop = RamseyPhase(pi,tPi,tPi.*[1/2],[720e-6],[0],X_col,Y_col,free_col,P_col,psi_col,1e6)
    pf_col[i] = pf_col_loop[1]
end
figure(2)
plot(Δ,abs.(pf_col .- 0.5).*2)
xlabel("Pk2Pk detuning")
ylabel("Ramsey at π phase")
title("Ramsey contrast with wait = 0.72ms")

# SIMULATE COLLECTIVE RAMSEY NO ECHO FOR 8 MOLECULES DETUNING ONLY TIME SCAN
num_molecules = 8;
rabi_frequencies =ones(8) #Fractional errors in pi time for each state
tWaits = range(0,5e-3,100)

b_coll = tensor([b for i=1:num_molecules]...)
Xs(i,rabRat) = embed(b_coll,i,Ω*σx*rabRat)
Ys(i,rabRat) = embed(b_coll,i,Ω*σy*rabRat)
frees(i,Δ) = embed(b_coll,i,Δ*P2)
P1s(i) = embed(b_coll,i,P1)
P_col = sum(P1s.(1:num_molecules))/num_molecules

X_col = sum(Xs.(1:num_molecules,rabi_frequencies))
Y_col = sum(Ys.(1:num_molecules,rabi_frequencies))
dets = 832.2*(range(-0.5,0.5,num_molecules)).*2*pi;
free_col = sum(frees.(1:num_molecules,dets))

pf_col = zeros(length(tWaits))

for i = 1:length(tWaits)
    psi_col = tensor([psi for i=1:num_molecules]...)
    pf_col_loop = RamseyPhase(pi,tPi,tPi.*[1/2],[tWaits[i]],[0],X_col,Y_col,free_col,P_col,psi_col,1e6)
    pf_col[i] = pf_col_loop[1]
end
figure(3)
plot(tWaits.*1e3,abs.(pf_col .- 0.5).*2)
xlabel("Ramsey time [ms]")
ylabel("Ramsey contrast at π phase")
=#
# SIMULATE COLLECTIVE XY8 AS A FUNCTION OF NGROUPS
#= maxGroups = 16;
tsXY = tPi.*vcat([1/2],ones(maxGroups*8))
tWaitsXY = ones(size(tsXY))*100*1e-6
phasesXY = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0],maxGroups))
phasesXX = vcat([0],repeat([0,0],maxGroups*4))
phasesXY16 = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0,pi,3*pi/2,pi,3*pi/2,3*pi/2,pi,3*pi/2,pi],Int(maxGroups/2)))
 =#
#= XYcontr = zeros(maxGroups)
XXcontr = zeros(maxGroups)
for i = 1:maxGroups
    res = RamseyPhase([pi],tPi,tsXY[1:8*i+1],tWaitsXY[1:8*i+1],phasesXY[1:8*i+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XYcontr[i] = res[end]
    resXX = RamseyPhase([pi],tPi,tsXY[1:8*i+1],tWaitsXY[1:8*i+1],phasesXX[1:8*i+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XXcontr[i] = resXX[end]
end

XY16contr = zeros(Int(maxGroups/2))
for i = 1:Int(maxGroups/2)
    res = RamseyPhase([pi],tPi,tsXY[1:8*i*2+1],tWaitsXY[1:8*i*2+1],phasesXY16[1:8*i*2+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XY16contr[i] = res[end]
end

figure(3)
plot(1:maxGroups,XXcontr,label = "X8")
plot(1:maxGroups,XYcontr,label = "XY8")
plot(2:2:maxGroups,XY16contr,label = "XY16")
xlabel("N groups")
ylabel("Ramsey contrast")
title("Avg pulse sequence contrasts with 1% Ω err and +- 500 Hz Δ err")
legend() =#

# Try adding some noise on specification of phases
#= phasesXYNoisy = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0],maxGroups))
phasesXYNoisy = phasesXYNoisy .+ pi/15*(rand(Int(length(phasesXYNoisy))) .- 0.5)
phasesXYNoisy[1] = 0;
print(phasesXYNoisy[1:10])

XYNoisycontr = zeros(maxGroups)
for i = 1:maxGroups
    res = RamseyPhase([pi],tPi,tsXY[1:8*i+1],tWaitsXY[1:8*i+1],phasesXYNoisy[1:8*i+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XYNoisycontr[i] = res[end]
end

figure(5)
plot(1:maxGroups,XYNoisycontr,label = "XY8 with ±π/30 phase noise")
xlabel("N groups")
ylabel("Ramsey contrast")
title("Avg pulse sequence contrasts with 1% Ω err and +- 500 Hz Δ err")
legend() =#

# Try an imperfect rotation to Y
#=phasesXYImperf = vcat([0],repeat([0,pi/2-pi/4,0,pi/2-pi/4,pi/2-pi/4,0,pi/2-pi/4,0],maxGroups))

XYImperfcontr = zeros(maxGroups)
for i = 1:maxGroups
    res = RamseyPhase([pi],tPi,tsXY[1:8*i+1],tWaitsXY[1:8*i+1],phasesXYImperf[1:8*i+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XYImperfcontr[i] = res[end]
end

figure(6)
plot(1:maxGroups,XYImperfcontr,label = "XY8 with Y axis off by π/4")
xlabel("N groups")
ylabel("Ramsey contrast")
title("Avg pulse sequence contrasts with 1% Ω err and +- 500 Hz Δ err")
legend()=#

# Plot XY8 state evolution
#= tsXY = tPi.*vcat([1/2],ones(80))
tWaitsXY = ones(size(tsXY))*100*1e-6
phasesXY = vcat([0],repeat([0,pi/2],40))
tout, psi_t = DDSeq(tsXY,tWaitsXY,phasesXY,Ω*σx,Ω*σy,Δ*P2,psi,1e6)
exp_val = expect(P1, psi_t)
exp_X = expect(σx, psi_t)
figure(6)
plot(tout,exp_val,label="N=0 popn")
plot(tout,exp_X,label="<ψ|σx|ψ>")
xlabel("Time (s)")
ylabel("Expectation value")
legend() =#
