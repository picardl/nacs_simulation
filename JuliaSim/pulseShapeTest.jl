using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using Profile

N = 2 #4
tPi = 1.6e-6/3; #Microwave pi pulse time at zero detuning
ΩBase = 2*pi*(1/(4*tPi));
dets = 2*pi*range(-0.905e6,-0.895e6,20);
#dets = 2*pi*range(-50e3,0,15);

exp_val_N1 =  Array{Any, 1}(undef, length(dets));
XRot, YRot, FreeEv, Ps, b = DynamicalDecoupling.genNLevelOperators(N, ΩBase, dets[1]);

#pulseShape = "square";
pulseShape = "truncGaussPulse";
#pulseShape = "truncGaussDiscrete";

for i=1:length(dets)
    print(".")

    params = Dict();
    params["frac"] = 0.3;
    params["len"] = 100e-6;

    psi = nlevelstate(b,1);

    amps = [1];
    times =[params["len"]];
    phases = [0];

    tTots,ψTots,Ps = DynamicalDecoupling.generalPulseSeq(psi,amps,times,phases,ΩBase,[dets[i]],N,pulseShape,params)

    exp_val_N1[i] = expect(Ps[2], ψTots[end])
end

figure(1)
plot(dets/(2*pi),exp_val_N1,label="N=1")
xlabel("Detuning (Hz)")
ylabel("Popn")
