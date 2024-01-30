using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using LsqFit
using Profile
using Random
using Dates
using DataFrames
using CSV

N = 2
tPi = 20.5e-6; #Microwave pi pulse time at zero detuning
inter_t = 1.32e-3;

h = 6.62607015e-34
hbar = h/(2*pi)
eps0 = 8.8541878128e-12
SR = 1e6;

NGauss = 100000;

dets = [500,600]

gaussStd = 2*pi*25;
gaussNoise = gaussStd*randn(NGauss) .+ 1
gaussFunc = LinearInterpolation(range(0,1,NGauss), gaussNoise)
 
figure(1)
plot(range(0,50e-3,5000),gaussFunc(range(0,50e-3,5000)))

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
    #2*pi*[dets[1]*fracIntens(t).+ gaussFunc(t)]; 
    2*pi*[dets[1]*(5*fracIntens(t)-4)]; 
end
Δt2 = (t) -> begin #Microwave detuning
    #2*pi*[500,-329.6e3,320e3]*fracIntens(t); 
    2*pi*[dets[2]]*(5*fracIntens(t)-4); 
end
Ωt = (t) -> begin
    #2*pi*(1/(4*tPi))*[1,2.805,0.58];
    #2*pi*(1/(4*tPi))*[1,0.2725,0.6585];
    2*pi*(1/(4*tPi))*[1];
end


XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperatorsTimeDep(N, Ωt, Δt1)
XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperatorsTimeDep(N, Ωt, Δt2)

XRotStat, YRotStat, FreeEvStat, PsStat, bStat = DynamicalDecoupling.genNLevelOperators(N, 2*pi*(1/(4*tPi))*[1], [2*pi*dets[1]])
FreeEv1 = (t,psi) -> Δt1(t)[1]*PsStat[2];

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

if true

#Single body spin echo
psi0 = nlevelstate(b1,1)
tWaits = range(10e-6,77e-3,10)
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = [0];#range(0,2*pi,15);

NTrials = 10;
mean_values = zeros(Float64,length(tWaits),1)
stderr_values = zeros(Float64,length(tWaits),1)
these_expect = zeros(Float64,1,NTrials)

println("Simulating single body coherence")
for i = 1:1:length(tWaits)
    print('.')
    tWaitsSpinEcho = tWaits[i]*DynamicalDecoupling.waitFracSpinEcho

    for j=1:NTrials
        tStart = rand(1)*0.2;
    
        local Xoffs = (t, psi) -> XRot1(t + tStart[1], psi)
        local Yoffs = (t,psi)-> YRot1(t + tStart[1], psi)
        local Freeoffs = (t,psi)-> FreeEv1(t + tStart[1], psi)

        pfN, ψfN= DynamicalDecoupling.RamseyPhaseTimeDep(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,Xoffs,Yoffs,Freeoffs,Ps1[1],psi0,SR)
        these_expect[:,j] = pfN
    end
    mean_values[i,:] = mean(these_expect, dims=2)
    stderr_values[i,:] = std(these_expect, dims=2)./sqrt(NTrials)
end
figure(2)
#plot(tWaits,expect(Ps1[1], ψfN),label="|0⟩")
errorbar(tWaits,mean_values[:,1],stderr_values[:,1])
legend()
xlabel("Ramsey time")
ylabel("Popn")
ylim([-0.01,1.01])
title("Single body contrast")
end

if true
#Spin echo time scan
psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
tWaits = vcat(range(10e-6,3*inter_t,20),range(12*inter_t,15*inter_t,20),range(60*inter_t,70*inter_t,60))
#tWaits = range(10e-6,15*inter_t,150)
tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
probePhases = 0

Ps = [P00,P01,P10,P11]

NTrials = 10;
mean_values = zeros(Float64,length(tWaits),4)
stderr_values = zeros(Float64,length(tWaits),4)
these_expect = zeros(Float64,4,NTrials)

println("Simulating interaction")
for i = 1:1:length(tWaits)
    print('.')
    for j=1:NTrials
        tStart = rand()*0.2;

        local Xoffs = (t, psi) -> X_col(t .+ tStart, psi)
        local Yoffs= (t,psi)-> Y_col(t .+ tStart, psi)
        local Freeoffs = (t,psi)-> free_col(t .+ tStart, psi)

        tWaitsSpinEcho = tWaits[i]*DynamicalDecoupling.waitFracSpinEcho
        pf, ψf = DynamicalDecoupling.RamseyPhaseTimeDep(probePhases,tPi,tsSpinEcho,tWaitsSpinEcho,DynamicalDecoupling.phasesSpinEcho,Xoffs,Yoffs,Freeoffs,P00,psi00,SR)
        these_expect[:,j] = [inner[1] for inner in real.(expect.(Ps, Ref(ψf)))]
    end
    mean_values[i,:] = mean(these_expect, dims=2)
    stderr_values[i,:] = std(these_expect, dims=2)./sqrt(NTrials)
end



figure(3)
errorbar(tWaits,mean_values[:,1],stderr_values[:,1],label="|00⟩")
errorbar(tWaits,mean_values[:,2],stderr_values[:,2],label="|01⟩")
errorbar(tWaits,mean_values[:,3],stderr_values[:,3],label="|10⟩")
errorbar(tWaits,mean_values[:,4],stderr_values[:,4],label="|11⟩")

legend()
xlabel("Spin echo wait time")
ylabel("Popn")
ylim([0,1])
title(current_datetime)

# Generate the CSV file name with the current date and time


# Save the data to a CSV file
data_dir = "data"
mkpath(data_dir)  # Create the 'data' directory if it doesn't exist
current_datetime = Dates.format(now(), "yyyymmdd_HHMMSS")
csv_filename = joinpath(data_dir, "data_$current_datetime.csv")

# Extract the data you want to save
expect_P00_values =mean_values[:,1]
expect_P01_values = mean_values[:,2]
expect_P10_values = mean_values[:,3]
expect_P11_values = mean_values[:,4]
data = hcat(tWaits, expect_P00_values, expect_P01_values,expect_P10_values,expect_P11_values)
data_df = DataFrame(tWaits=data[:, 1], Expect_P00=data[:, 2], Expect_P01=data[:, 3],Expect_P10=data[:, 4],Expect_P11=data[:, 5])

# Write the data to the CSV file
CSV.write(csv_filename, data_df, header=["tWaits", "Expect_P00",  "Expect_P01", "Expect_P10","Expect_P11"])

end