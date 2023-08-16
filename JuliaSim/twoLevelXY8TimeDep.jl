using QuantumOptics
using PyPlot
using Interpolations
using Statistics
using LsqFit

function DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR,tOffs)
    #Does time evolution for a generic dynamical decoupling pulse sequence
    #ts: (N,1) array of rotation times 
    #tWaits: (N,1) array of free evolution times after each pulse in ts
    #phases: (N,1) array of phases of each rotation
    #XRot: Operator representing rotation about X axis (rabi freq should already by factor in)
    #YRot: Operator representing rotation about Y axis (rabi freq should already by factor in)
    #freeEv: Operator representing free evolution between pulses, which should also capture detuning of states during pulse
    #ψ0: initial state Vector
    #SR: sample rate for time points
    ψ = ψ0
    T = 1/SR
    tEnd = tOffs

    for i = 1:length(ts)
        tPulse = tEnd:T:ts[i]+tEnd
        H_t = (t,psi) -> begin
            freeEv(t,psi) + cos(phases[i])*XRot + sin(phases[i])*YRot
        end
        tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, ψ, H_t)

        tEnd = tout[end]
        ψ = ψ_t[end]

        if tWaits[i] > 0
            tPulse = tEnd:T:tWaits[i]+tEnd
            tout, ψ_t = timeevolution.schroedinger_dynamic(tPulse, ψ, freeEv)

            tEnd = tout[end]
            ψ = ψ_t[end]
        end
    end

    return tEnd, ψ
end

function RamseyPhase(probePhases,tPi,ts,tWaits,phases,XRot,YRot,freeEv,P1,ψ0,SR,tOffs)
    #Return final |0⟩ popn for final pi/2 phases in probePhases
    tEnd, ψ = DDSeq(ts,tWaits,phases,XRot,YRot,freeEv,ψ0,SR,tOffs)
    exp_val = expect(P1, ψ)
    
    pf = zeros(size(probePhases));
    for i = 1:length(probePhases)
        H_t = (t,psi) -> begin
            freeEv(t,psi) + cos(probePhases[i])*XRot + sin(probePhases[i])*YRot
        end
        toutf, psif = timeevolution.schroedinger_dynamic([0,tPi/2] .+ tEnd, ψ,H_t)
        pf[i] = real(expect(P1, psif[end]))
    end
    return pf
end

csv_file_path = "C:/nacs_simulation/JuliaSim/dcNoise250ms.CSV"
#csv_file_path = "C:/nacs_simulation/JuliaSim/8Traps5s.CSV"
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

N = 2
b = NLevelBasis(N)
Δ = 2*pi*(0); #Microwave detuning
tPi = 32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi))

psi = nlevelstate(b,1)

maxGroups = 32
tsXY = tPi.*vcat([1/2],ones(maxGroups*8))
tWaitsXY = ones(size(tsXY))*4e-3
phasesXY = vcat([0],repeat([0,pi/2,0,pi/2,pi/2,0,pi/2,0],maxGroups))
 
#Sigma x and y operators for the two levels
σx = (transition(b,1,2) + dagger(transition(b,1,2)))
σy = (-im*transition(b,1,2) + im*dagger(transition(b,1,2)))

#Projection operators for each state
P1 = tensor(nlevelstate(b,1), dagger(nlevelstate(b,1)))
P2 =  tensor(nlevelstate(b,2), dagger(nlevelstate(b,2)))

#general rotation Hamiltonian
genRot(Ωarg,ϕ) = Ωarg*(cos(ϕ)*σx + sin(ϕ)*σy)

# SIMULATE COLLECTIVE RAMSEY SPIN ECHO
num_molecules = 8;
#rabi_frequencies =range(0.995,1.005,num_molecules) #Fractional errors in pi time for each state
rabi_frequencies = [49.6,49.5,48.5,49.9,49.9,49.9,49.9,50.1]./mean([49.6,49.5,48.5,49.9,49.9,49.9,49.9,50.1])

b_coll = tensor([b for i=1:num_molecules]...)
Xs(i,rabRat) = embed(b_coll,i,Ω*σx*rabRat)
Ys(i,rabRat) = embed(b_coll,i,Ω*σy*rabRat)
frees(i,Δ) = embed(b_coll,i,Δ*P2)
P1s(i) = embed(b_coll,i,P1)
P_col = sum(P1s.(1:num_molecules))/num_molecules

X_col = sum(Xs.(1:num_molecules,rabi_frequencies))
Y_col = sum(Ys.(1:num_molecules,rabi_frequencies))
#dets = 832.2*(range(-0.5,0.5,num_molecules)).*2*pi;
dets = ([-78.2,14.0,108.6,186.2,152.3,339.7,339.7,508.0]).*2*pi;
#dets = rand(8).*2*pi.*781;

free_col = sum(frees.(1:num_molecules,dets))

free_col_t = (t,psi) -> begin
    free_col * 1#(fracIntens(t))
end

psi_col = tensor([psi for i=1:num_molecules]...)
nGroups = [1,4]
tWaits = nGroups*16*4e-3
cts = zeros(length(tWaits))
ctsErrs = zeros(length(tWaits))
NShots = 20
NPhases = 20
pf_shots = zeros(NShots,NPhases)
@time begin
    for j = 1:length(nGroups)
        alltOffs = rand(1,NShots)
        print('.') 
        for i = 1:NShots
            pf_col = RamseyPhase(range(0,2*pi,NPhases),tPi,tsXY[1:8*nGroups[j]+1],tWaitsXY[1:8*nGroups[j]+1],phasesXY[1:8*nGroups[j]+1],X_col,Y_col,free_col_t,P_col,psi_col,1e6,alltOffs[i])
            pf_shots[i,:] = pf_col;
#=              figure(j)
            plot(range(0,2*pi,30),pf_col)  =#
        end
        pf_avg = mean(pf_shots,dims = 1)
        pf_se = std(pf_shots,dims = 1)./sqrt(NShots)
        sMax =  maximum(pf_avg)
        sMaxErr = pf_se[argmax(pf_avg)]
        sMin =  minimum(pf_avg)
        sMinErr =  pf_se[argmin(pf_avg)]
        cts[j] = sMax -sMin
        ctsErrs[j] = sMaxErr + sMinErr
        figure(1)
        plot(range(0,2*pi,NPhases),pf_avg[1,:])  
        xlabel("Ramsey phase")
        ylabel("N=0 popn")
    end
end
 ftfun(x,p) = p[1].*exp.(-x.^2.0/abs(p[2]).^2.0) .+ p[3]
guess = [1.0,0.6,0.6]
ftResult = curve_fit(ftfun,tWaits,cts,1 ./ctsErrs,guess)
fit_params = ftResult.param
#fit_errs = stderror(ftResult) 

figure(77)
errorbar(tWaits,cts,ctsErrs)
xfit = range(0,tWaits[end],1000)
plot(xfit,ftfun(xfit,fit_params))
print("Decay time: $(fit_params[2])")
xlabel("Spin echo time [s]")
ylabel("Ramsey coherence")
ylim([0,1])
#= pf_col = zeros(length(tWaits))

for i = 1:length(tWaits)
    psi_col = tensor([psi for i=1:num_molecules]...)
    pf_col_loop = RamseyPhase(pi,tPi,tsSpinEcho,tWaitsSpinEcho,phasesSpinEcho,X_col,Y_col,free_col_t,P_col,psi_col,1e6)
    pf_col[i] = pf_col_loop[1]
end
figure(3)
plot(tWaits.*1e3,abs.(pf_col .- 0.5).*2)
xlabel("Ramsey time [ms]")
ylabel("Ramsey contrast at π phase") =#

