using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using LsqFit
using Profile
using CSV
using Dates
using DataFrames

N = 4
tPi = 20.5e-6; #Microwave pi pulse time at zero detuning
inter_t = 1.32e-3;

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12

#csv_file_path = "C:/nilab-projects/nacs_simulation/JuliaSim/dcNoise250ms.CSV"
csv_file_path = "C:/projects/nacs_simulation/JuliaSim/dcNoise250ms.CSV"
lines = readlines(csv_file_path)
lines = map(strip, lines)
data = map(line -> split(line, ','), lines)
# Convert the 4th and 5th columns to Float64, leave other columns as strings
data = map(row -> [row[1], row[2], row[3], parse(Float64, row[4]), parse(Float64, row[5]), row[6:end]], data)
# Extract the 4th column
data_t = map(row -> row[4], data)
# Extract the 5th column
data_v = map(row -> row[5], data)
data_t = data_t[100:end] .- data_t[100]
data_v = data_v[100:end]

fracIntens = LinearInterpolation(data_t, data_v./mean(data_v))

Δt1 = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    2*pi*[500,-261.17e3,261.17e3]*fracIntens(t); 
end
Δt2 = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    2*pi*[600,-261.17e3,261.17e3]*fracIntens(t); 
end
Ωt = (t) -> begin
    #2*pi*(1/(4*tPi))*[1,2.805,0.58];
    #2*pi*(1/(4*tPi))*[1,0.2725,0.6585];
    2*pi*(1/(4*tPi))*[1,1,1];
end


XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperatorsTimeDep(N, Ωt, Δt1)
XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperatorsTimeDep(N, Ωt, Δt2)

b_col = b1 ⊗ b2

X_col = (t,psi) -> embed(b_col,1,XRot1(t,psi)) + embed(b_col,2,XRot2(t,psi))
Y_col = (t,psi) -> embed(b_col,1,YRot1(t,psi)) + embed(b_col,2,YRot2(t,psi))

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

free_col = (t,psi) -> embed(b_col,1,FreeEv1(t,psi)) + embed(b_col,2,FreeEv2(t,psi)) + H_int


#Spin echo time scan
psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
tWaits = [10e-6,inter_t];#vcat([10e-6,inter_t],range(11*inter_t,13*inter_t,2))
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = 0
pfN = Array{Any, 1}(undef, length(tWaits));
ψfN = Array{Any, 1}(undef, length(tWaits));

for i = 1:1:length(tWaits)
    print('.')
    tWaitsSpinEcho = tWaits[i]*DynamicalDecoupling.waitFracSpinEcho
    pfN[i], ψfN[i] = DynamicalDecoupling.RamseyPhaseTimeDep(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,X_col,Y_col,free_col,P00,psi00,1e6)
end

figure(3)
plot(tWaits,expect(P00, ψfN),label="N=0")
plot(tWaits,expect(P11, ψfN),label="N=1,0")

legend()
xlabel("Spin echo wait time")
ylabel("Popn")
ylim([0,1])

# Save the data to a CSV file
data_dir = "data"
mkpath(data_dir)  # Create the 'data' directory if it doesn't exist

# Generate the CSV file name with the current date and time
current_datetime = Dates.format(now(), "yyyymmdd_HHMMSS")
csv_filename = joinpath(data_dir, "data_$current_datetime.csv")

# Extract the data you want to save
expect_P00_values = [real(value[1]) for value in expect(P00, ψfN)]
expect_P11_values = [real(value[1])  for value in expect(P11, ψfN)]
data = hcat(tWaits, expect_P00_values, expect_P11_values)
data_df = DataFrame(tWaits=data[:, 1], Expect_P00=data[:, 2], Expect_P11=data[:, 3])

# Write the data to the CSV file
CSV.write(csv_filename, data_df, header=["tWaits", "Expect_P00", "Expect_P11"])
