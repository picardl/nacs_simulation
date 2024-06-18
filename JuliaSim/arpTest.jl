using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using Profile

N = 5;
tPi_HF = 310e-6;
tPi = tPi_HF/193; #Microwave pi pulse time at zero detuning

Δe = 0.720e6;
Δg = 1.89e6;

Ωs = 2*pi*(1/(4*tPi))*[1,1,0.86];

Δs = [-Δg,0,-Δe,0,Δe]



if false #ARP

    df = 25e3;
rampRate = 0.1e9;

params = Dict();
params["frac"] = 0.3;
params["len"] = df/rampRate;
shapeName = "square";
#shapeName = "truncGaussPulse";
#shapeName = "truncGaussDiscrete";



b = NLevelBasis(N);

psi0 = nlevelstate(b,2);

Ωt = DynamicalDecoupling.loadPulseShape(shapeName,1 .*Ωs,params);
ΩHFt = DynamicalDecoupling.loadPulseShape(shapeName,Ωs./193,params);

    Δt = (t) -> Δs .+  (-df/2 .+ df*t/params["len"]).*[0,0,1,1,1];

XRot  = (t,psi) -> begin
    return sum([Ω * (transition(b, 1, i+2) + dagger(transition(b, 1, i+2))) for (i, Ω) in enumerate(Ωt(t))]) + sum([Ω *(transition(b, 2, i+2) + dagger(transition(b, 2, i+2))) for (i, Ω) in enumerate(ΩHFt(t))])
end
YRot  = (t,psi) -> begin
    return sum([Ω * (-im*transition(b, 1, i+1) + im*dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωt(t))]) + sum([Ω * (-im*transition(b, 2, i+2) + im*dagger(transition(b, 2, i+2))) for (i, Ω) in enumerate(ΩHFt(t))])
end

Ps = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:N];

FreeEv  = (t,psi) -> begin
    return sum([Δ * Ps[i] for (i, Δ) in enumerate(Δt(t))])
end

tTots = []
ψTots = []
    

    amps = [1];
    times =range(1e-6,params["len"],10);
    phases = [0];
    tEnd = 0;

    for i = 1:length(times)
        print(".")
        tPulse = [0,times[i]];
        

        H_t = (t,psi) -> begin
            FreeEv(t,psi) + amps[1]*cos(phases[1])*XRot(t,psi)
        end

        tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, psi0,  H_t);
        global tTots = vcat(tTots,tout[end] .+ tEnd);
        global ψTots = vcat(ψTots,ψ_t[end]);

        global tEnd = tTots[end];
        ψ = ψ_t[end];
    end

    exp_val_N1 = expect(Ps[4],ψTots)

figure(1)
plot(times,expect(Ps[1],ψTots),label="|0⟩")
plot(times,expect(Ps[2],ψTots),label="|1⟩")
plot(times,expect(Ps[3],ψTots),label="|N=1,lower⟩")
plot(times,expect(Ps[4],ψTots),label="|e⟩")
plot(times,expect(Ps[5],ψTots),label="|N=1,upper⟩")
xlabel("Time")
ylabel("Popn")
legend()

end

if true #Rabi det


params = Dict();
params["frac"] = 0.3;
params["len"] = tPi_HF;
shapeName = "square";
#shapeName = "truncGaussPulse";
#shapeName = "truncGaussDiscrete";



b = NLevelBasis(N);

psi0 = nlevelstate(b,2);

Ωt = DynamicalDecoupling.loadPulseShape(shapeName,1 .*Ωs,params);
ΩHFt = DynamicalDecoupling.loadPulseShape(shapeName,Ωs./193,params);


XRot  = (t,psi) -> begin
    return sum([Ω * (transition(b, 1, i+2) + dagger(transition(b, 1, i+2))) for (i, Ω) in enumerate(Ωt(t))]) + sum([Ω *(transition(b, 2, i+2) + dagger(transition(b, 2, i+2))) for (i, Ω) in enumerate(ΩHFt(t))])
end
YRot  = (t,psi) -> begin
    return sum([Ω * (-im*transition(b, 1, i+1) + im*dagger(transition(b, 1, i+1))) for (i, Ω) in enumerate(Ωt(t))]) + sum([Ω * (-im*transition(b, 2, i+2) + im*dagger(transition(b, 2, i+2))) for (i, Ω) in enumerate(ΩHFt(t))])
end

Ps = [tensor(nlevelstate(b, i), dagger(nlevelstate(b, i))) for i in 1:N];



tTots = []
ψTots = []
    

    amps = [1];
    dets = 2*pi*range(-60e3,-50e3,15);
    phases = [0];
    tEnd = 0;

    for i = 1:length(dets)
        Δt = (t) -> Δs .+  dets[i].*[0,0,1,1,1];
        FreeEv  = (t,psi) -> begin
            return sum([Δ * Ps[i] for (i, Δ) in enumerate(Δt(t))])
        end


        print(".")
        tPulse = [0,params["len"]];
        

        H_t = (t,psi) -> begin
            FreeEv(t,psi) + amps[1]*cos(phases[1])*XRot(t,psi)
        end

        tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, psi0,  H_t);
        global tTots = vcat(tTots,tout[end] .+ tEnd);
        global ψTots = vcat(ψTots,ψ_t[end]);

        global tEnd = tTots[end];
        ψ = ψ_t[end];
    end

    exp_val_N1 = expect(Ps[4],ψTots)

figure(1)
plot(dets./(2*pi),expect(Ps[1],ψTots),label="|0⟩")
plot(dets./(2*pi),expect(Ps[2],ψTots),label="|1⟩")
plot(dets./(2*pi),expect(Ps[3],ψTots),label="|N=1,lower⟩")
plot(dets./(2*pi),expect(Ps[4],ψTots),label="|e⟩")
plot(dets./(2*pi),expect(Ps[5],ψTots),label="|N=1,upper⟩")
xlabel("Detuning / 2Pi (Hz)")
ylabel("Popn")
legend()

end