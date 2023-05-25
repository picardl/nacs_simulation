using QuantumOptics
using PyPlot

N = 2
b = NLevelBasis(N)
Δ = 2*pi*(0); #Microwave detuning
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

micro = Ω*(transition(b,1,2) + dagger(transition(b,1,2))) #Microwave coupling

tspan = range(0,8e-3,1000000)

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

#Evaluate multi-particle version by defining a collective Hamiltonian
# Initialize an empty list to store the individual molecule states
num_molecules = 8;
rabi_frequencies = range(0.995,1.005,8) #Fractional errors in pi time for each state


b_coll = tensor([b for i=1:num_molecules]...)
Hs(i,rabRat) = embed(b_coll,i,Δ*P2 + micro*rabRat) #Microwave drive Hamiltonian
P1s(i) = embed(b_coll,i,P1)

H_col = sum(Hs.(1:num_molecules,rabi_frequencies))
P_col = sum(P1s.(1:num_molecules))/num_molecules

psiCol = tensor([psi for i=1:num_molecules]...)
Tout, psiColT = timeevolution.schroedinger(tspan, psiCol, H_col)
expValCol = expect(P_col, psiColT)
figure(3)
plot(Tout,expValCol)
xlabel("Time")
ylabel("Average |0⟩ popn")