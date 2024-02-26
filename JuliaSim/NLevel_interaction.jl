using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12

N = 4
tPi = 20e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))*[1,0.5,0.5]
#Ω = 2*pi*(1/(4*tPi))*[1]

inter_t = 2.5e-3;

dets1 = [-25,700e3,700e3]
dets2 = [25,700e3,700e3]

XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets1)
XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets2)
#XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets)
#XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets)
b_col = b1 ⊗ b2

#Xs(i) = embed(b2,i,XRot)
#Ys(i) = embed(b2,i,YRot)

X_col = embed(b_col,1,XRot1) + embed(b_col,2,XRot2)
Y_col = embed(b_col,1,YRot1) + embed(b_col,2,YRot2)
free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2)

fudge = 1/(2*pi*(J/(4*hbar)))/2

θ = 0/180*pi;
R = 2e-6;
d = 4.6*3.33564e-30
J = (d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2)

fudge = abs(1/(J/(2*hbar))/2*pi)/inter_t

H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)))
P00 = Ps1[1] ⊗ Ps2[1];
P01= Ps1[1] ⊗ Ps2[2];
P10= Ps1[2] ⊗ Ps2[1];
P11= Ps1[2] ⊗ Ps2[2];

P0 = Ps1[1] ⊗ identityoperator(b2);

# Vanilla interaction starting in 01
#=
psi01 = nlevelstate(b1,1) ⊗ nlevelstate(b2,2)
tout, psi_t = timeevolution.schroedinger(range(0,0.05,1000), psi01, H_int + free_col)
exp_val = expect(P01, psi_t)
plot(tout,exp_val,label="|00⟩")
plot(tout,expect(P0, psi_t),label="|0⟩")
xlabel("Time (s)")
ylabel("Popn in |01>")
title("Interaction starting in 01")=#

# Rabi
#=
figure(2)
psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
tout, psi_t = timeevolution.schroedinger(range(0,100e-6,1000), psi00, Y_col + free_col)
exp_val = expect(P00, psi_t)
plot(tout,exp_val,label="|00⟩")
plot(tout,expect(P0, psi_t),label="|0⟩")
xlabel("Time (s)")
ylabel("Popn in |01>")
title("Rabi")
legend()
=#


#Ramsey spin echo interaction starting in 00

psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
#tWaits = vcat(range(1e-6,5e-3,60),range(20e-3,25e-3,60))
tWaits = range(1e-6,5e-3,60)
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = 0
pfN = Array{Any, 1}(undef, length(tWaits));
ψfN = Array{Any, 1}(undef, length(tWaits));

for i = 1:1:length(tWaits)
    tWaitsSpinEcho = tWaits[i]*DynamicalDecoupling.waitFracSpinEcho
    pfN[i], ψfN[i] = DynamicalDecoupling.RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,X_col,Y_col, free_col + H_int,P00,psi00,1e6)
end

figure(3)
plot(tWaits,real(expect(P00, ψfN)),label="|00⟩")
plot(tWaits,real(expect(P01, ψfN)),label="|01⟩")
plot(tWaits,real(expect(P10, ψfN)),label="|10⟩")
plot(tWaits,real(expect(P11, ψfN)),label="|11⟩")
legend()
xlabel("Spin echo wait time")
ylabel("Popn")
ylim([0,1.01]);


# Ramsey spin echo phases
#=
psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
tWaits = 100e-6;
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = range(0,2*pi,30)
pfN = Array{Any, 1}(undef, length(probePhases));
ψfN = Array{Any, 1}(undef, length(probePhases));

tWaitsSpinEcho = tWaits*DynamicalDecoupling.waitFracSpinEcho

pfN, ψfN = DynamicalDecoupling.RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,X_col,Y_col, free_col,P00,psi00,1e6)

figure(4)
plot(probePhases,real(expect(P00,ψfN)),label="|00⟩")
plot(probePhases,real(expect(P01, ψfN)),label="|01⟩")
plot(probePhases,real(expect(P0, ψfN)),label="|0⟩")
legend()
xlabel("Probe phase")
ylabel("Popn")
ylim([0,1.01]);
=#