module DynamicalDecoupling

using QuantumOptics
using PyPlot

export DDSeq, RamseyPhase, genNLevelOperators,genXY8, tFracSpinEcho, waitFracSpinEcho, phasesSpinEcho

function DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    """
    DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    Perform time evolution for a generic dynamical decoupling pulse sequence.

    Arguments:
    - ts::Vector: (N,1) array of rotation times 
    - tWaits::Vector: (N,1) array of free evolution times after each pulse in ts
    - phases::Vector: (N,1) array of phases of each rotation
    - XRot: Operator representing rotation about X axis (rabi freq should already by factor in)
    - YRot: Operator representing rotation about Y axis (rabi freq should already by factor in)
    - freeEv: Operator representing free evolution between pulses, which should also capture detuning of states during pulse
    - ψ0::Vector: initial state Vector
    - SR::Float64 sample rate for time points

    Returns:
    - tTots::Vector: List of times
    - ψTots::Vector: List of state vectors at each time point
    """

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
    """
    RamseyPhase(probePhases, tPi, ts, tWaits, phases, XRot, YRot, freeEv, P1, ψ0, SR)

    Perform Ramsey phase experiment and return final |0⟩ population for given probe phases.
    
    Arguments:
    - probePhases::Vector: List of probe phases for the Ramsey experiment.
    - tPi::Float64: Duration of the π pulse.
    - ts::Vector: List of rotation times for the DD sequence.
    - tWaits::Vector: List of free evolution times after each pulse in `ts`.
    - phases::Vector: List of phases of each rotation in `ts`.
    - XRot: Operator representing rotation about X axis (rabi freq should already be factored in).
    - YRot: Operator representing rotation about Y axis (rabi freq should already be factored in).
    - freeEv: Operator representing free evolution between pulses, capturing detuning of states during pulse.
    - P1: Projection operator for the |0⟩ state.
    - ψ0::Vector: Initial state vector.
    - SR::Float64: Sample rate for time points.
    
    Returns:
    - pf::Vector: List of final |0⟩ populations corresponding to each probe phase.
    - ψf::Vector: Vector of final state vectors after tine evoltuion
    """

    tout, psi_t = DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    exp_val = expect(P1, psi_t)
    
    pf = zeros(size(probePhases));
    ψf = Array{Any, 1}(undef, length(probePhases))
    for i = 1:length(probePhases)
        toutf, psif = timeevolution.schroedinger([0,tPi/2], psi_t[end], freeEv + cos(probePhases[i])*XRot + sin(probePhases[i])*YRot)
        pf[i] = real(expect(P1, psif[end]))
        ψf[i] = psif[end]
    end
    return pf, ψf
end

function CollectiveRamseyPhase(probePhases,tPi,ts,tWaits,phases,XRot,YRot,freeEv,P1,ψ0,SR)
    """
    CollectiveRamseyPhase(probePhases, tPi, ts, tWaits, phases, XRot, YRot, freeEv, P1, ψ0, SR)

    Perform a set of Ramsey phase experiments and return average final |0⟩ population for given probe phases.
    
    Arguments:
    - probePhases::Vector: List of probe phases for the Ramsey experiment.
    - tPi::Float64/Vector: Duration of the π pulse. Either single value or Vector.
    - ts::Vector: List of rotation times for the DD sequence.
    - tWaits::Vector: List of free evolution times after each pulse in `ts`.
    - phases::Vector: List of phases of each rotation in `ts`.
    - XRot/Vector: Operator representing rotation about X axis (rabi freq should already be factored in). Either single value or Vector.
    - YRot/Vector: Operator representing rotation about Y axis (rabi freq should already be factored in). Either single value or Vector.
    - freeEv/Vector: Operator representing free evolution between pulses, capturing detuning of states during pulse. Either single value or Vector.
    - P1: Projection operator for the |0⟩ state.
    - ψ0::Vector: Initial state vector.
    - SR::Float64: Sample rate for time points.
    
    Returns:
    - pf::Vector: List of final |0⟩ populations corresponding to each probe phase.
    - pf::Vector: Vector of final state vectors after tine evoltuion
    """

    allRes = RamseyPhase.(Ref(probePhases),tPi,Ref(ts),Ref(tWaits),Ref(phases),XRot,YRot,freeEv,Ref(P1),Ref(ψ0),SR)

    numExperiments = length(allRes)  # Number of experiments

    # Loop to average over all the |0⟩ state populations for the different possible inputs
    avgPf = sum([res[1] for res in allRes]) / numExperiments
    
    return avgPf, allRes
end

function genXY8(tPi,tau,nGroups)
    """
    genXY8(nGroups)
    Generate XY8 sequence parameters for a given number of groups

    Arguments:
    - tPi::Float64: Pi pulse time
    - tPi::tau: tau time between pi pulses
    - nGroups::Integer: number of groups of XY8 pulses

    
    Returns:
    - tsXY::Vector: List of pulse times
    - tWaitsXY::Vector: List of wait times
    - phasesXY::Vector: List of pulses phases
    """

    tsXY = tPi.*vcat([1/2],ones(nGroups*8))
    tWaitsXY = ones(size(tsXY))*tau
    phasesXY = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0],nGroups))
    return tsXY,tWaitsXY,phasesXY
end

function genNLevelOperators(N, Ωs, Δs)
    """
    genNLevelOperators(N, Ωs, Δs)
    Generate rotation, free evolution, and projection operators for N level basis.
    Assumed that rotation is between level 1 and each of the other levels, with no rotations between higher levels.

    Arguments:
    - N::Integer: Number of levels
    - Ωs::Vector: Vector of length N-1, Rabi frequency of rotation from level 1 to each other level
    - Δs::Vector: Vector of length N-1, detunings of each level in microwave rotating frame

    Returns:
    - XRot::Operator: Rotation operator about X axis for all levels
    - YRot::Operator: Rotation operator about Y axis for all levels
    - FreeEv::Operator: Sum of free evolution operators for all levels
    - Ps::Vector: Array of projection operators onto each of the N levels
    """
    b = NLevelBasis(N)
    
    XRot = sum([Ω * (transition(b, 1, i+1) + dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωs)])
    YRot = sum([Ω * (-im*transition(b, 1, i+1) + im*dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωs)])
    
    Ps = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:N]
    FreeEv = sum([Δ * Ps[i+1] for (i, Δ) in enumerate(Δs)])

    return XRot, YRot, FreeEv, Ps
end

tFracSpinEcho = [1/2,1]
waitFracSpinEcho = [1,1]
phasesSpinEcho = [0,0]


end;