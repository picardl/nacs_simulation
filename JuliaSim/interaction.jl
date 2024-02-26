using QuantumOptics
using PyPlot

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12

N = 2
b = NLevelBasis(N)
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

dets = 2*pi*[-100,100]

θ = 0/180*pi;
R = 2e-6;
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
P01 = P0 ⊗ P1;
P10 = P1 ⊗ P0;
P11 = P1 ⊗ P1;
Xs(i) = embed(b2,i,Ω*σx)
Ys(i) = embed(b2,i,Ω*σy)
frees(i,Δ) = embed(b2,i,Δ*P1)
P0s(i) = embed(b2,i,P0)

X_col = sum(Xs.(1:2))
Y_col = sum(Ys.(1:2))
free_col = sum(frees.(1:2,dets))

# Vanilla interaction starting in 01
#psi01 = nlevelstate(b,1) ⊗ nlevelstate(b,2)
#tout, psi_t = timeevolution.schroedinger(range(0,0.02,1000), psi01, H_int + free_col)
#exp_val = expect(P0s(1), psi_t)
#plot(tout,exp_val)


# Ramsey spin echo interaction starting in 00
psi00 = nlevelstate(b,1) ⊗ nlevelstate(b,1)
tWaits = range(1e-6,20e-3,100)
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = 0
pfN = Array{Any, 1}(undef, length(tWaits));
ψfN = Array{Any, 1}(undef, length(tWaits));

for i = 1:1:length(tWaits)
    tWaitsSpinEcho = tWaits[i]*DynamicalDecoupling.waitFracSpinEcho
    pfN[i], ψfN[i] = DynamicalDecoupling.RamseyPhase(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,X_col,Y_col, H_int + free_col,P00,psi00,1e6)
end

figure(3)
plot(tWaits,expect(P00, ψfN),label="|00⟩")
plot(tWaits,expect(P01, ψfN),label="|01⟩")
plot(tWaits,expect(P10, ψfN),label="|10⟩")
plot(tWaits,expect(P11, ψfN),label="|11⟩")
legend()
xlabel("Spin echo wait time")
ylabel("Popn")
ylim([0,1])