using QuantumOptics
using PyPlot
using Random

N = 2
b = NLevelBasis(N)
Δ = 2*pi*(0); #Microwave detuning
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

micro = Ω*(transition(b,1,2) + dagger(transition(b,1,2))) #Microwave coupling

tspan = range(0,8e-3,100000)

#Diagonal terms
P1 = tensor(nlevelstate(b,1), dagger(nlevelstate(b,1)))
P2 =  tensor(nlevelstate(b,2), dagger(nlevelstate(b,2)))

H = Δ*P2 + micro #Microwave drive Hamiltonian

psi = nlevelstate(b,1)

#Evolve single particle Hamiltonian
tout, psi_t = timeevolution.schroedinger(tspan, psi, H)
exp_val = expect(P1, psi_t)
figure(1)
plot(tout,exp_val)
xlabel("Time (s)")
ylabel("Average |0⟩ popn")

#Evaluate multi-particle version by defining a collective Hamiltonian
# Initialize an empty list to store the individual molecule states
num_molecules = 8;
#rabi_frequencies = range(0.995,1.005,8) #Fractional errors in pi time for each state
rabi_frequencies = ones(8) #Fractional errors in pi time for each state

b_coll = tensor([b for i=1:num_molecules]...)
Hs(i,rabRat,Δ) = embed(b_coll,i,Δ*P2 + micro*rabRat) #Microwave drive Hamiltonian
P1s(i) = embed(b_coll,i,P1)

H_col = sum(Hs.(1:num_molecules,rabi_frequencies,Δ))
P_col = sum(P1s.(1:num_molecules))/num_molecules

psiCol = tensor([psi for i=1:num_molecules]...)
Tout, psiColT = timeevolution.schroedinger(tspan, psiCol, H_col)
expValCol = expect(P_col, psiColT)
figure(3)
plot(Tout,expValCol,label="Δ = 0")
xlabel("Time")
ylabel("Average |0⟩ popn")

#Check effect of adding detuning as well
#dets = range(-2e3,2e3,8)*2*pi
dets = 500*(rand(8).-0.5).*2*pi
H_col_det = sum(Hs.(1:num_molecules,rabi_frequencies,dets))

Tout, psiColT = timeevolution.schroedinger(tspan, psiCol, H_col_det)
expValCol = expect(P_col, psiColT)
figure(3)
plot(Tout,expValCol,label="Δ = -500 to 500 Hz")
xlabel("Time")
ylabel("Average |0⟩ popn")
legend()
