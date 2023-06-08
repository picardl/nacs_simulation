using QuantumOptics
using PyPlot

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12

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
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

dets = [-100,100]

θ = 35/180*pi;
R = 2.6e-6;
d = 4.6*3.33564e-30
J = (d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2)

b2 = b ⊗ b

#Sigma x and y operators for the two levels
σx = (transition(b,1,2) + dagger(transition(b,1,2)))
σy = (-im*transition(b,1,2) + im*dagger(transition(b,1,2)))

H_int = J/hbar/2*(transition(b,1,2) ⊗ transition(b,2,1) + dagger(transition(b,1,2) ⊗ transition(b,2,1)))

#Projection operators for each state
P0 = tensor(nlevelstate(b,1), dagger(nlevelstate(b,1)))
P1 =  tensor(nlevelstate(b,2), dagger(nlevelstate(b,2)))

P00 = P0 ⊗ P0;

Xs(i) = embed(b2,i,Ω*σx)
Ys(i) = embed(b2,i,Ω*σy)
frees(i,Δ) = embed(b2,i,Δ*P1)
P0s(i) = embed(b2,i,P0)

X_col = sum(Xs.(1:2))
Y_col = sum(Ys.(1:2))
free_col = sum(frees.(1:2,dets))

# Vanilla interaction starting in 01
psi01 = nlevelstate(b,1) ⊗ nlevelstate(b,2)
tout, psi_t = timeevolution.schroedinger(range(0,0.1,1000), psi01, H_int + free_col)
exp_val = expect(P0s(1), psi_t)
plot(tout,exp_val)