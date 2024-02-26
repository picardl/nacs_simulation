module DynamicalDecoupling
using Interpolations
using QuantumOptics
using PyPlot
using SpecialFunctions
#using StochasticDiffEq

export DDSeq, RamseyPhase,CollectiveRamseyPhase
export DDSeqTimeDep,RamseyPhaseTimeDep,CollectiveRamseyPhaseTimeDep
export genNLevelOperators,genNLevelOperatorsTimeDep,genXY8, tFracSpinEcho, waitFracSpinEcho, phasesSpinEcho

function generalPulseSeq(ψ0,amps,times,phases,Omega,Deltas,N,pulseShape,params)

    tTots = []
    ψTots = []
    ψ = ψ0;
    tEnd = 0;

    _, _, _, Ps, _ = genNLevelOperators(N, 1, 1)

    for i = 1:length(times)
        Δt = (t) -> Deltas; 
        params["len"] = times[i]
        Ωt = loadPulseShape(pulseShape,Omega,params);
        XRot, YRot, FreeEv,_, _ = genNLevelOperatorsTimeDep(N, Ωt, Δt);
        tPulse = [0,times[i]];

        H_t = (t,psi) -> begin
            FreeEv(t,psi) + amps[i]*cos(phases[i])*XRot(t,psi) + amps[i]*sin(phases[i])*YRot(t ,psi)
        end

        tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, ψ,  H_t);
        tTots = vcat(tTots,tout[end] .+ tEnd);
        ψTots = vcat(ψTots,ψ_t[end]);

        tEnd = tTots[end];
        ψ = ψ_t[end];
    end

    return tTots,ψTots,Ps

end

function generalPulseSeq_stochastic(ψ0,amps,times,phases,Omega,Deltas,H_noise,N,pulseShape,params)

    tTots = []
    ψTots = []
    ψ = ψ0;
    tEnd = 0;
    dt = 10e-9;

    _, _, _, Ps, _ = genNLevelOperators(N, 1, 1)

    for i = 1:length(times)
        Δt = (t) -> Deltas; 
        params["len"] = times[i]
        Ωt = loadPulseShape(pulseShape,Omega,params);
        XRot, YRot, FreeEv,_, _ = genNLevelOperatorsTimeDep(N,  Ωt, Δt);
        tPulse = [0,times[i]];

        H_t = (t,psi) -> begin
            return FreeEv(t,psi) + amps[i]*cos(phases[i])*XRot(t,psi) + amps[i]*sin(phases[i])*YRot(t ,psi), [H_noise], [dagger(H_noise)]
        end
        # H_noise_t = (t,psi) -> begin
        #     [H_noise]
        # end
        H = FreeEv(0,0) + amps[i]*cos(phases[i])*XRot(0,0) + amps[i]*sin(phases[i])*YRot(0,0)

        tout, ψ_t = stochastic.schroedinger_dynamic(tPulse, ψ,  H_t,H_noise_t;alg=StochasticDiffEq.RKMilGeneral(interpretation=:Stratonovich),dt=dt,abstol = 1e-3);
        #tout, ψ_t = timeevolution.master_dynamic(tPulse, ψ,  H_t,rates = [1]);
        #tout, ψ_t = timeevolution.master(tPulse, ψ,  H,[H_noise];rates = [1]);

        tTots = vcat(tTots,tout[end] .+ tEnd);
        ψTots = vcat(ψTots,ψ_t[end]);

        tEnd = tTots[end];
        ψ = ψ_t[end];
    end

    return tTots,ψTots,Ps

end

function generalPulseSeq_master(ψ0,amps,times,phases,Omega,Delta,H_noise,N,params)

    tTots = []
    rhoTots = []
    rho = ψ0⊗dagger(ψ0);
    tEnd = 0;
    dt = 10e-9;

    XRot, YRot, FreeEv,Ps, _ = genNLevelOperators(N,  Omega, Delta);

    for i = 1:length(times)
        Δt = (t) -> Deltas; 
        params["len"] = times[i]
        
        tPulse = [0,times[i]];

        H = FreeEv + amps[i]*cos(phases[i])*XRot + amps[i]*sin(phases[i])*YRot

        tout, rho_t = timeevolution.master(tPulse, rho,  H,[H_noise];rates = [1]);

        tTots = vcat(tTots,tout[end] .+ tEnd);
        rhoTots = vcat(rhoTots,rho_t[end]);

        tEnd = tTots[end];
        rho = rho_t[end];
    end

    return tTots,rhoTots,Ps

end

function loadPulseShape(name,Omega,params)
    if name == "square"
        Ωt = (t) -> Omega
    elseif name == "truncGaussPulse"
        frac = params["frac"];
        len = params["len"];
        tau = frac*len;
        intPulse = sqrt(pi)*tau*erf(len/(2*tau));
        intRect = len*Omega;
        ampMax = intRect/intPulse;
        Ωt = (t) -> ampMax*exp(-(t-len/2).^2/tau.^2);
    elseif name == "truncGaussDiscrete"
        frac = params["frac"];
        len = params["len"];
        tau = frac*len;
        intPulse = sqrt(pi)*tau*erf(len/(2*tau));
        intRect = len*Omega;
        ampMax = intRect/intPulse;
        Ωt = (t) -> ampMax*exp(-(ceil(t/2e-6)*2e-6-len/2).^2/tau.^2);
    end
    return Ωt
end

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
        ψ = ψ_t[end]

        if tWaits[i] > 0
            tPulse = 0:T:tWaits[i]
            tout, ψ_t = timeevolution.schroedinger(tPulse, ψ, freeEv)
            tTots = vcat(tTots,tout .+ tEnd)
            ψTots = vcat(ψTots,ψ_t)

            tEnd = tTots[end]
            ψ = ψ_t[end]
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

function DDSeqTimeDep(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
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

    T = 1/SR
    
    tTots = Array{Float64}(undef, Integer(sum(floor.(ts/T).+1) + sum(floor.(tWaits/(10*T)).+1)))
    ψTots = Array{Any}(undef, Integer(sum(floor.(ts/T).+1) + sum(floor.(tWaits/(10*T)).+1)))
    ψ = ψ0
    tEnd = 0
    indEnd = 0
    for i = 1:length(ts)
        tPulse = 0:T:ts[i]

        H_t = (t,psi) -> begin
            freeEv(t + tEnd,psi) + cos(phases[i])*XRot(t + tEnd,psi) + sin(phases[i])*YRot(t + tEnd,psi)
        end

        tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, ψ, H_t)
        tTots[indEnd+1:indEnd+length(tout)] = tout .+ tEnd
        ψTots[indEnd+1:indEnd+length(tout)] = ψ_t

        indEnd = indEnd+length(tout)
        tEnd = tTots[indEnd]
        ψ = ψTots[indEnd]

        if tWaits[i] > 0
            tPulse = 0:T*10:tWaits[i]
            tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, ψ, (t,psi) -> freeEv(t + tEnd,0))
            tTots[indEnd+1:indEnd+length(tout)] = tout .+ tEnd
            ψTots[indEnd+1:indEnd+length(tout)] = ψ_t
            indEnd = indEnd+length(tout)
            tEnd = tTots[indEnd]
            ψ = ψTots[indEnd]
        end
    end
    return tTots, ψTots
end

function RamseyPhaseTimeDep(probePhases,tPi,ts,tWaits,phases,XRot,YRot,freeEv,P1,ψ0,SR)
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

    tout, psi_t = DDSeqTimeDep(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR)
    exp_val = expect(P1, psi_t[end])
    pf = zeros(size(probePhases));
    ψf = Array{Any, 1}(undef, length(probePhases))
    for i = 1:length(probePhases)
        H_t = (t,psi) -> begin
            freeEv(t + tout[end],psi) + cos(probePhases[i])*XRot(t+ tout[end],psi) + sin(probePhases[i])*YRot(t+ tout[end],psi) 
        end
        toutf, psif = timeevolution.schroedinger_dynamic([0,tPi/2], psi_t[end], H_t)
        pf[i] = real(expect(P1, psif[end]))
        ψf[i] = psif[end]
    end
    return pf, ψf
end

function CollectiveRamseyPhaseTimeDep(probePhases,tPi,ts,tWaits,phases,XRot,YRot,freeEv,P1,ψ0,SR)
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

    allRes = RamseyPhaseTimeDep.(Ref(probePhases),tPi,Ref(ts),Ref(tWaits),Ref(phases),XRot,YRot,freeEv,Ref(P1),Ref(ψ0),SR)

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

    tsXY = tPi.*vcat([1/2],ones(nGroups*8));
    tWaitsXY = ones(size(tsXY))*tau*2;
    tWaitsXY[1] = tau;
    tWaitsXY[length(tWaitsXY)] = tau;
    phasesXY = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0],nGroups));
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
    - b::NLevel: N level basis object
    """
    b = NLevelBasis(N)
    
    XRot = sum([Ω * (transition(b, 1, i+1) + dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωs)])
    YRot = sum([Ω * (-im*transition(b, 1, i+1) + im*dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωs)])
    
    Ps = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:N]
    FreeEv = sum([Δ * Ps[i+1] for (i, Δ) in enumerate(Δs)])

    return XRot, YRot, FreeEv, Ps, b
end

function genNLevelOperatorsTimeDep(N, Ωt, Δt)
    """
    genNLevelOperators(N, Ωs, Δs)
    Generate rotation, free evolution, and projection operators for N level basis.
    Assumed that rotation is between level 1 and each of the other levels, with no rotations between higher levels.

    Arguments:
    - N::Integer: Number of levels
    - Ωs::Function: Function of t, returning vector of length N-1, Rabi frequency of rotation from level 1 to each other level
    - Δs::Function: Function of t, returning vector of length N-1, detunings of each level in microwave rotating frame

    Returns:
    - XRot::Function: Rotation operator about X axis for all levels as a function of t and psi (unused argument)
    - YRot::Function: Rotation operator about Y axis for all levels as a function of t and psi (unused argument)
    - FreeEv::Function: Sum of free evolution operators for all levels as a function of t and psi (unused argument)
    - Ps::Vector: Array of projection operators onto each of the N levels
    - b::NLevel: N level basis object
    """
    b = NLevelBasis(N);
    
    X = LazySum([0.0 for i in 1:N-1],[(transition(b, 1, i+1) + dagger(transition(b, 1, i+1))) for i in 1:N-1]);
    Y = LazySum([0.0 for  i in 1:N-1],[(-im*transition(b, 1, i+1) + im*dagger(transition(b, 1, i+1))) for  i in 1:N-1]);

    XRot  = (t,psi) -> begin
        X.factors = ([Ω for (i, Ω) in enumerate(Ωt(t))]);
        return X
    end
    YRot  = (t,psi) -> begin
        Y.factors = ([Ω for (i, Ω) in enumerate(Ωt(t))]);
        return Y
    end
    
    Ps = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:N];

    Ev = LazySum([0.0 for i in 1:N-1],[Ps[i+1] for  i in 1:N-1]);

    FreeEv  = (t,psi) -> begin
        Ev.factors =[Δ for (i, Δ) in enumerate(Δt(t))];
        return Ev
    end

    return XRot, YRot,FreeEv, Ps, b
end

tFracSpinEcho = [1/2,1];
waitFracSpinEcho = [1,1];
phasesSpinEcho = [0,0];


end;