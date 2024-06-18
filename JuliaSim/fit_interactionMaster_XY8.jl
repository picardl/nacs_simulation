
using QuantumOptics
using PyPlot
using Optimization
using OptimizationBBO
using OptimizationOptimJL
using CSV
using DataFrames
#using Zygote
using ForwardDiff
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Random
#using Plots

function masterSurvivals(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 2
    tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];
    tau = params["tau"];
    nGroups = params["nGroups"];

    Ω = 2*pi*(1/(4*tPi))*[1];

    inter_t = abs(x[1]);

    dets1 = [-abs(x[2])/2];
    dets2 = [abs(x[2])/2];
    gammaCOM = abs(x[3]);
    gammaREL = abs(x[4]);


    J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = abs(x[5]);

    scalingFactor = abs(x[6])

    XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets1);
    XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets2);
    b_col = b1 ⊗ b2;

    X_col = embed(b_col,1,XRot1) + embed(b_col,2,XRot2);
    Y_col = embed(b_col,1,YRot1) + embed(b_col,2,YRot2);
    free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2);

    H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    P00 = Ps1[1] ⊗ Ps2[1];
    P01= Ps1[1] ⊗ Ps2[2];
    P10= Ps1[2] ⊗ Ps2[1];
    P11= Ps1[2] ⊗ Ps2[2];

    int_col = free_col + H_int

    P0 = Ps1[1] ⊗ identityoperator(b2);

    #H_noise_1B = sum([sqrt(gamma) * Ps1[i+1] for (i,gamma) in enumerate(gammas1)])⊗sum([sqrt(gamma2) * Ps2[j+1] for (j,gamma2) in enumerate(gammas2)]);
    H_noise_DeltaCOM = sqrt(gammaCOM)*(-1*P00 + P11);
    H_noise_DeltaREL = sqrt(gammaREL/2)*(-1*P01 + P10);
    H_noise_2B = sqrt(gammaJ)*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    H_noise = [H_noise_DeltaCOM, H_noise_DeltaREL,H_noise_2B];

    psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)
    probePhases = 0
    rhof = Array{Any, 1}(undef, length(tData));
    
    for i = 1:1:length(tData)
        amps,ts,phases = DynamicalDecoupling.genXY8(tPi,tau,nGroups[i])
        phases[end] = probePhases
        tTots,rhoTots = DynamicalDecoupling.generalTwoBodyPulseSeq_master(psi00,amps,ts,phases,X_col,Y_col,int_col,H_noise,params);
        rhof[i] = rhoTots[end];
    end

    fitData = zeros(4,length(rhof));
    fitData[1,:] = scalingFactor*real(expect(P00, rhof));
    fitData[2,:] = scalingFactor*real(expect(P01, rhof));
    fitData[3,:] = scalingFactor*real(expect(P10, rhof));
    fitData[4,:] = scalingFactor*real(expect(P11, rhof));

    return fitData
end

#=function masterResid(x,params)
    tData = params["tData"];
    yData = params["yData"];
    errLower = params["errLower"];
    errUpper = params["errUpper"];

    nGroupsInit = params["nGroups"];

    newParams = copy(params)

    resids = zeros(1,8)
    for i = 1:8
        randInd = rand(1:length(nGroupsInit))
        newParams["nGroups"] = [nGroupsInit[randInd]]
        fitData = masterSurvivals(x,tData[randInd],newParams);

        weights = 1.0./(errLower[randInd].^2.0 + errUpper[randInd].^2.0 .+ 1e-9);
        resids[i] = sum(((fitData - yData[:,randInd]).^2.0).*weights);
    end
    print(sum(resids))
    print("\n")
    return sum(resids)

end=#

function masterResid(x,params)
    tData = params["tData"];
    yData = params["yData"];
    errLower = params["errLower"];
    errUpper = params["errUpper"];

    fitData = masterSurvivals(x,tData,params);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((fitData - yData).^2.0).*weights;
    print(sum(resids))
    print("\n")
    return sum(resids)
end


function callback_function(opt_values,f_val)
    # Extracting the iteration count and function value
    global iter += 1;
    #f_val = opt_values.f_val

    # Plotting the iteration and function value
    if iter%10 == 0
        figure(1)
        semilogy(iter, f_val; marker="o", color="red")
        xlabel("Iteration")
        ylabel("Function Value")
        title("Optimization Progress")
        draw()  # Update the plot in real-time
        pause(0.01)
    end
        # Plotting the iteration and function value
        if iter%10 == 0
            figure(2)
            for i = 1:length(opt_values)
                subplot(length(opt_values),1,i)
                plot(iter, opt_values[i]; marker="o", color="red")
                xlabel("Iteration")
                ylabel("Value")
                draw()  # Update the plot in real-time
            end
        end
    return false
end

#=
function load_and_extract_data(filename)
    # Load CSV file into a DataFrame
    df = CSV.File(filename) |> DataFrame

    # Separate data into two sets based on scanIdx
    data_scan1 = filter(row -> row.scanIdx == 1, df)
    data_scan2 = filter(row -> row.scanIdx == 2, df)

    # Extract vectors for each set
    groupIdx_scan1 = data_scan1.groupIdx
    axisValue_scan1 = data_scan1.axisValue
    mean_scan1 = data_scan1.mean
    errorLower_scan1 = data_scan1.errorLower
    errorUpper_scan1 = data_scan1.errorUpper

    groupIdx_scan2 = data_scan2.groupIdx
    axisValue_scan2 = data_scan2.axisValue
    mean_scan2 = data_scan2.mean
    errorLower_scan2 = data_scan2.errorLower
    errorUpper_scan2 = data_scan2.errorUpper

    return (axisValue_scan1, mean_scan1, errorLower_scan1, errorUpper_scan1,
            axisValue_scan2, mean_scan2, errorLower_scan2, errorUpper_scan2)
end
=#

function load_and_extract_data(filename)
    # Load CSV file into a DataFrame
    df = CSV.File(filename) |> DataFrame

    # Extract vectors for each set
    axisValue_scan = df.axisValue
    mean_scan = df.mean
    errorLower_scan = df.errorLower
    errorUpper_scan = df.errorUpper
    nGroups_scan = df.nGroups

    return (axisValue_scan, mean_scan, errorLower_scan, errorUpper_scan,nGroups_scan)
end

function stirapNorm(stirap_contrast,stirap_contrastErrLower,stirap_contrastErrUpper,this_surv,this_errLower,this_errUpper)
    numerator = (this_surv .-  stirap_contrast[1]);
    stirap_diff2B = (stirap_contrast[2] - stirap_contrast[1])
    stirap_diff2BErrLower = sqrt(stirap_contrastErrLower[1].^2 + stirap_contrastErrLower[2].^2);
    stirap_diff2BErrUpper = sqrt(stirap_contrastErrUpper[1].^2 + stirap_contrastErrUpper[2].^2);
    this_surv = (this_surv .-  stirap_contrast[1])/stirap_diff2B;
    this_errNumUpper = sqrt.(this_errUpper.^2 .+ stirap_contrastErrUpper[1].^2);
    this_errNumLower = sqrt.(this_errLower.^2 .+ stirap_contrastErrLower[1].^2);
    this_errLower = this_surv.*sqrt.((this_errNumLower./numerator).^2 .+ (stirap_diff2BErrLower./stirap_diff2B).^2);
    this_errUpper = this_surv.*sqrt.((this_errNumUpper./numerator).^2 .+ (stirap_diff2BErrUpper./stirap_diff2B).^2);
    this_errLower[isnan.(this_errLower)] .= this_errNumLower[isnan.(this_errLower)]./stirap_diff2B; 
    this_errUpper[isnan.(this_errUpper)] .= this_errNumUpper[isnan.(this_errUpper)]./stirap_diff2B; 
    return (this_surv,this_errLower,this_errUpper)
end

global iter = 0;

dStamp = "20240402"
tStamp = "200923"

#dataPath = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/20240225_115611"
dataPath = "C:/projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

(testt, survival00, errLower00, errUpper00,nGroups) = load_and_extract_data(dataPath*"/00_data"*dStamp*"_"*tStamp*".csv")
(_,survival01, errLower01, errUpper01,_) = load_and_extract_data(dataPath*"/01_data"*dStamp*"_"*tStamp*".csv")
(_,survival10, errLower10, errUpper10,_) = load_and_extract_data(dataPath*"/10_data"*dStamp*"_"*tStamp*".csv")
(_,survival11, errLower11, errUpper11,_) = load_and_extract_data(dataPath*"/11_data"*dStamp*"_"*tStamp*".csv")

#=
survival00, errLower00, errUpper00 = stirapNorm(stirap00,stirapErrLower00,stirapErrUpper00,survival00,errLower00,errUpper00);
survival01, errLower01, errUpper01 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival01, errLower01, errUpper01)
survival10, errLower10, errUpper10 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival10, errLower10, errUpper10)
survival11, errLower11, errUpper11 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival11, errLower11, errUpper11)
=#

allSurvival = hcat(survival00, survival01, survival10, survival11)'
allErrLower = hcat(errLower00, errLower01, errLower10, errLower11)'
allErrUpper = hcat(errUpper00, errUpper01, errUpper10, errUpper11)'

params = Dict();
params["tPi"] = 20e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
params["θ"] = 0/180*pi;
params["R"] = 2e-6;
params["d"] = 4.6*3.33564e-30
params["inter_t"] = 2e-3; #Interaction pi time
Delta = 2*pi*50; #Site-by-site detuning in 2*pi*Hz
params["tau"] = testt[1]/nGroups[1]/16
params["nGroups"] = nGroups


guess = [1.4e-3,Delta,2*pi*10,2*pi*10,2*pi*10,0.8]; #interaction time, detuning, COM detuning noise, Relative detuning noise, interaction noise, normalization factor (scales all points)

#testt = [1e-6,1.5e-3,3e-3,7.5e-3,9e-3]
#testData = [0.06 0.015 0.045 0.05 0.02;0 0.01 0.002 0.002 0.005;0 0.01 0.002 0.002 0.005;0 0.025 0.005 0.005 0.02;]./0.06;
#testErrLower = [0.01 0.005 0.01 0.01 0.01;0 0.005 0.002 0.002 0.005;0 0.005 0.002 0.002 0.005;0 0.005 0.002 0.002 0.005]./0.06;
#testErrUpper = [0.01 0.005 0.01 0.01 0.01;0.002 0.005 0.002 0.002 0.005;0.002 0.005 0.002 0.002 0.005;0.002 0.005 0.002 0.002 0.005]./0.06;

params["tData"] = testt;
params["yData"] = allSurvival;
params["errLower"] = allErrLower;
params["errUpper"] = allErrUpper;
#prob = OptimizationProblem(masterResid, guess,params, lb = [0.5e-3,0,0,0,0], ub = [3e-3,2*pi*1e3,2*pi*1e3,2*pi*1e3,2*pi*1e3])
prob = OptimizationProblem(masterResid, guess,params)

#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); callback = callback_function, local_reltol = 1e-3,local_abstol = 1e-3,maxiters = 10000)
#sol = solve(prob, ParticleSwarm(); callback = callback_function,x_tol = 1e-5, f_tol = 1e-3);
#sol = solve(prob, SAMIN(); callback = callback_function, x_tol = 1e-6, f_tol = 1e-3)
sol = solve(prob, NelderMead(); callback = callback_function, x_tol = 1e-6, f_tol = 1e-3,maxiters = 5000)

print(sol)

#resid = masterResid(guess,params)
#print(resid);

nPlot = range(1,maximum(nGroups));
tPlot = 16*params["tau"]*nPlot;
params["nGroups"] = nPlot;
out = masterSurvivals(sol,tPlot,params);

figure(3);
plot(tPlot,out[1,:],label="|00⟩",color="C0");
errorbar(testt,allSurvival[1,:],yerr=hcat(allErrLower[1,:],allErrUpper[1,:])',color="C0",linestyle="none",marker="o")
plot(tPlot,out[2,:],label="|01⟩",color="C1");
errorbar(testt,allSurvival[2,:],yerr=hcat(allErrLower[2,:],allErrUpper[2,:])',color="C1",linestyle="none",marker="o")
plot(tPlot,out[3,:],label="|10⟩",color="C2");
errorbar(testt,allSurvival[3,:],yerr=hcat(allErrLower[3,:],allErrUpper[3,:])',color="C2",linestyle="none",marker="o")
plot(tPlot,out[4,:],label="|11⟩",color="C3");
errorbar(testt,allSurvival[4,:],yerr=hcat(allErrLower[4,:],allErrUpper[4,:])',color="C3",linestyle="none",marker="o")

legend();
xlabel("Spin echo wait time");
ylabel("Popn");
ylim([0,1.01]);