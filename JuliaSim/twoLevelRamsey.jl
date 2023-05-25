using QuantumOptics
using PyPlot

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

N = 2
b = NLevelBasis(N)
Δ = 2*pi*(0e3); #Microwave detuning
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

tspan = range(0,100e-6,10000)

tpi2 =range(0,tPi/2,1000)

#Sigma x and y operators for the two levels
σx = (transition(b,1,2) + dagger(transition(b,1,2)))
σy = (-im*transition(b,1,2) + im*dagger(transition(b,1,2)))

#Projection operators for each state
P1 = tensor(nlevelstate(b,1), dagger(nlevelstate(b,1)))
P2 =  tensor(nlevelstate(b,2), dagger(nlevelstate(b,2)))

#general rotation Hamiltonian
genRot(Ωarg,ϕ) = Ωarg*(cos(ϕ)*σx + sin(ϕ)*σy)

# SIMULATE SPIN ECHO FOR A SINGLE MOLECULE
tsSpinEcho = tPi.*[1/2,1]
tWaitsSpinEcho = [100,100]*1e-6
phasesSpinEcho = [0,pi]

probePhases = range(0,2*pi,30)

pf = RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,Ω*σx,Ω*σy,Δ*P2,P1,psi,1e6)
figure(1)
plot(probePhases,pf)
xlabel("Ramsey phase")
ylabel("N=0 popn")
title("Single molecule spin Echo Ramsey")

# SIMULATE COLLECTIVE SPIN ECHO FOR 8 MOLECULES
num_molecules = 8;
rabi_frequencies = range(0.995,1.005,8) #Fractional errors in pi time for each state
dets = range(-500,500,8)

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
title("Avg Spinecho Ramsey with 1% Ω err and +- 500 Hz Δ err")

# SIMULATE COLLECTIVE XY8 AS A FUNCTION OF NGROUPS

tsXY = tPi.*vcat([1/2],ones(80))
tWaitsXY = ones(size(tsXY))*100*1e-6
phasesXY = vcat([0],repeat([0,pi/2],40))
phasesXX = vcat([0],repeat([0,0],40))

XYcontr = zeros(10)

for i = 1:10
    res = RamseyPhase([pi],tPi,tsXY[1:8*i+1],tWaitsXY[1:8*i+1],phasesXY[1:8*i+1],X_col,Y_col,free_col,P_col,psi_col,1e6)
    XYcontr[i] = res[end]
end

figure(3)
plot(1:10,XYcontr)
xlabel("N XY8 groups")
ylabel("Ramsey contrast")
title("Avg XY8 contrast with 1% Ω err and +- 500 Hz Δ err")