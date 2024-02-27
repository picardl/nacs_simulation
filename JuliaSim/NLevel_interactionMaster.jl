#Simulate molecule-molecule interactions incorporating N levels (N>=2)
#using master equation with two terms to model decoherence

using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12

N = 2#4
tPi = 20e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
#Ω = 2*pi*(1/(4*tPi))*[1,0.5,0.5]
Ω = 2*pi*(1/(4*tPi))*[1]

inter_t = 2.5e-3;

#dets1 = 2*pi* [-25,700e3,700e3]
#dets2 = 2*pi*[25,700e3,700e3]
dets1 = 2*pi* [-200]
dets2 = 2*pi*[200]
gammas1 = 0.01*abs.(dets1);
gammas2 = 0.01*abs.(dets2);

θ = 0/180*pi;
R = 2e-6;
d = 4.6*3.33564e-30
J = (d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2)
fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
gammaJ = 0*abs(fudge*J/hbar/2);

XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets1)
XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets2)
b_col = b1 ⊗ b2

#Xs(i) = embed(b2,i,XRot)
#Ys(i) = embed(b2,i,YRot)

X_col = embed(b_col,1,XRot1) + embed(b_col,2,XRot2)
Y_col = embed(b_col,1,YRot1) + embed(b_col,2,YRot2)
free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2)


H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)))
P00 = Ps1[1] ⊗ Ps2[1];
P01= Ps1[1] ⊗ Ps2[2];
P10= Ps1[2] ⊗ Ps2[1];
P11= Ps1[2] ⊗ Ps2[2];

int_col = free_col + H_int

P0 = Ps1[1] ⊗ identityoperator(b2);

H_noise_1B = sum([sqrt(gamma) * Ps1[i+1] for (i,gamma) in enumerate(gammas1)])⊗sum([sqrt(gamma2) * Ps2[j+1] for (j,gamma2) in enumerate(gammas2)]);
H_noise_2B = sqrt(gammaJ)*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
H_noise = H_noise_1B + H_noise_2B;


if true
    psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)

    amps = [1,0,1,0,1];
    phases::Vector{Float64} = [0,0,0,0,0];

    params = Dict();

    tWaits = range(1e-6,20e-3,60)
    tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
    probePhases = 0
    rhof = Array{Any, 1}(undef, length(tWaits));
    
    for i = 1:1:length(tWaits)
        print(".")
        ts = [tPi/2,tWaits[i]/2,tPi,tWaits[i]/2,tPi/2]
        tTots,rhoTots = DynamicalDecoupling.generalTwoBodyPulseSeq_master(psi00,amps,ts,phases,X_col,Y_col,int_col,H_noise,params);
        rhof[i] = rhoTots[end];
    end

    figure(3)
    plot(tWaits,real(expect(P00, rhof)),label="|00⟩")
    plot(tWaits,real(expect(P01, rhof)),label="|01⟩")
    plot(tWaits,real(expect(P10, rhof)),label="|10⟩")
    plot(tWaits,real(expect(P11, rhof)),label="|11⟩")
    legend()
    xlabel("Spin echo wait time")
    ylabel("Popn")
    ylim([0,1.01]);
end
