
using QuantumOptics
using PyPlot
using Optimization
using OptimizationBBO
using OptimizationOptimJL
using CSV
using DataFrames
using FiniteDiff
using Distributions
using HypothesisTests
#using Zygote
#using ForwardDiff
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
#using Plots

function masterSurvivals(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 2
    #=tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];=#
    tPi = params[1];##Microwave pi pulse time at zero detuning
    θ = params[2]; #θ
    R  = params[3]; #Spacing R
    d = params[4]; #dipole moment d

    Ω = 2*pi*(1/(4*tPi))*[1];

    inter_t = abs(x[1]);#Interaction pi time

    dets1 = [x[2]/2];
    dets2 = [x[2]/2];
    gammaCOM = abs(x[3]);
    gammaREL = abs(x[4]);


    J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = abs(x[5]);

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

    amps = [1,0,1,0,1];
    phases::Vector{Float64} = [0,0,pi/2,0,0];

    tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
    probePhases = 0

    rhof = Array{Any, 1}(undef, length(tData));
    
    for i = 1:1:length(tData)
        ts = [tPi/2,tData[i]/2,tPi,tData[i]/2,tPi/2]
        tTots,rhoTots = DynamicalDecoupling.generalTwoBodyPulseSeq_master(psi00,amps,ts,phases,X_col,Y_col,int_col,H_noise,params);
        rhof[i] = rhoTots[end];
    end

    fitData = zeros(4,length(rhof));
    fitData[1,:] = real(expect(P00, rhof));
    fitData[2,:] = real(expect(P01, rhof));
    fitData[3,:] = real(expect(P10, rhof));
    fitData[4,:] = real(expect(P11, rhof));

    return fitData
end

function masterResid(x,params)
    tData = params[6];
    yData = params[7];
    errLower = params[8];
    errUpper = params[9];

    fitData = masterSurvivals(x,tData,params);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((fitData - yData).^2.0).*weights;
    print(sum(resids))
    print("\n")
    return sum(resids)

end

function masterResidGlob(x)
    tData = globparams[6];
    yData = globparams[7];
    errLower = globparams[8];
    errUpper = globparams[9];

    fitData = masterSurvivals(x,tData,globparams);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((fitData - yData).^2.0).*weights;
    return reshape(resids,(1,116))

end

#=
function masterResid(x,params)
    #=tData = params["tData"];
    yData = params["yData"];
    errLower = params["errLower"];
    errUpper = params["errUpper"];=#
    tData = params[6];
    yData = params[7];
    errLower = params[8];
    errUpper = params[9];

    randInd = rand(1:length(tData),12)
    tSampled = [tData[i] for i in randInd]
    errLowerSampled = [errLower[i] for i in randInd]
    errUpperSampled = [errUpper[i] for i in randInd]

    yDataSampled = zeros(size(yData,1),12)
    for i = 1:length(randInd)
        yDataSampled[:,i] = yData[:,randInd[i]]
    end

    fitData = masterSurvivals(x,tSampled,params);
    weights = 1.0./(errLowerSampled.^2.0 + errUpperSampled.^2.0 .+ 1e-9);
    resids = sum(((fitData - yDataSampled).^2.0).*transpose(weights));

    print(sum(resids))
    print("\n")
    return sum(resids)

end=#

function callback_function(opt_values,f_val)
    # Extracting the iteration count and function value
    global iter += 1;
    #f_val = opt_values.f_val

    # Plotting the iteration and function value
    #=if iter%20 == 0
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
        end =#
    return false
end


function load_and_extract_data(filename)
    # Load CSV file into a DataFrame
    df = CSV.File(filename) |> DataFrame

    # Extract vectors for each set
    axisValue_scan = df.axisValue
    mean_scan = df.mean
    errorLower_scan = df.errorLower
    errorUpper_scan = df.errorUpper

    return (axisValue_scan, mean_scan, errorLower_scan, errorUpper_scan)
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

function genBootData(survival00,survival01,survival10,survival11,trials)
    nPoints = length(survival00)
    allSurvival = zeros(4,nPoints)
    allErrLower = zeros(4,nPoints)
    allErrUpper = zeros(4,nPoints)

    samples = zeros(4,nPoints)

    dists = Vector{Any}(undef,4)
    survs = [survival00,survival01,survival10,survival11]
    
    for i = 1:nPoints
        for j = 1:4
            dists[j] = Distributions.Binomial(trials[i],survs[j][i]);
            samples[j,i] = rand(dists[j]);
            conf = confint(HypothesisTests.BinomialTest(samples[j,i],trials[i]);level=0.682);

            allSurvival[j,i] = samples[j,i]/trials[i];
            allErrLower[j,i] = allSurvival[j,i] - conf[1];
            allErrUpper[j,i] =  conf[2] - allSurvival[j,i];
        end
    end
    return(allSurvival,allErrLower,allErrUpper)
end

global iter = 0;

#dStamp = "20240419"
#tStamp = "190010"
#dStamp = "20240501" #2.6 um data
#tStamp = "10534"
dStamp = "20240503" #2.3 um data
tStamp = "132113"
dataPath = "C:/projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

xScale = 1e-3;

(testt, survival00, errLower00, errUpper00) = load_and_extract_data(dataPath*"/00_data"*dStamp*"_"*tStamp*".csv")
(_,survival01, errLower01, errUpper01) = load_and_extract_data(dataPath*"/01_data"*dStamp*"_"*tStamp*".csv")
(_,survival10, errLower10, errUpper10) = load_and_extract_data(dataPath*"/10_data"*dStamp*"_"*tStamp*".csv")
(_,survival11, errLower11, errUpper11) = load_and_extract_data(dataPath*"/11_data"*dStamp*"_"*tStamp*".csv")

testt = testt*xScale;

df = CSV.File(dataPath*"/rawTrials_data"*dStamp*"_"*tStamp*".csv") |> DataFrame
s00 = df."00";
s01 = df."01";
s10 = df."10";
s11 = df."11";
trials = df.trials;

NBoots = 300;
solBoots = zeros(5,NBoots)

for bootInd = 1:NBoots

    (allSurvival,allErrLower,allErrUpper) =  genBootData(survival00,survival01,survival10,survival11,trials)

    Delta = 2*pi*50; #Site-by-site detuning in 2*pi*Hz

    params = (Vector{Any}(undef,9))
    params[1] = 13.4e-6;##Microwave pi pulse time at zero detuning
    params[2] = 0/180*pi; #θ
    params[3] = 2e-6; #Spacing R
    params[4] = 4.6*3.33564e-30 #dipole moment d
    params[5] = 2e-3; #Interaction pi time
    Delta = 2*pi*50; #Site-by-site detuning in 2*pi*Hz

    tPlot = range(1e-6,maximum(testt),1000);

    guess = [2e-3, 1000, 2*pi*2, 2*pi*2,2*pi*20]; #interaction time, detuning, COM detuning noise, Relative detuning noise, interaction noise

    params[6] = testt;
    params[7] = allSurvival;
    params[8] = allErrLower;
    params[9] = allErrUpper;
    

    paramTuple = Tuple(x for x in params)

    prob = OptimizationProblem(masterResid, guess,paramTuple)
    sol = solve(prob, NelderMead(lower = [1e-3,-2*pi*1000,0,0,0],upper=[4e-3,2*pi*1000,2*pi*100,2*pi*100,2*pi*100]); callback = callback_function, x_tol = 1e-7, f_tol = 1e-7,maxiters = 5000)
    print(sol)

    solBoots[:,bootInd] = sol;


    #out = masterSurvivals(sol,tPlot,paramTuple);

    #edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
    #faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

    #iSorted = sortperm(testt)
    #xPlotScale = 1e3;

    #=
    rc("font",family="sans",size = 14)
    figure(3);
    plot(tPlot*xPlotScale,out[1,:],label="|00⟩",color=edgeColors[1],lineWidth = 2);
    #errorbar(testt,allSurvival[1,:],yerr=hcat(allErrLower[1,:],allErrUpper[1,:])',color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
    plot(testt*xPlotScale,allSurvival[1,:],color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",markerSize = 7)
    fill_between(testt[iSorted]*xPlotScale,allSurvival[1,iSorted] - allErrLower[1,iSorted],allSurvival[1,iSorted] + allErrUpper[1,iSorted],color = faceColors[1],alpha = 0.3)

    plot(tPlot*xPlotScale,out[2,:],label="|01⟩",color=edgeColors[2],lineWidth = 2);
    #errorbar(testt,allSurvival[2,:],yerr=hcat(allErrLower[2,:],allErrUpper[2,:])',color = edgeColors[2],mec=edgeColors[2],mfc=faceColors[2],linestyle="none",marker="o",capsize = 3)
    plot(testt*xPlotScale,allSurvival[2,:],color = edgeColors[2],mec=edgeColors[2],mfc=faceColors[2],linestyle="none",marker="o",markerSize = 7)
    fill_between(testt[iSorted]*xPlotScale,allSurvival[2,iSorted] - allErrLower[2,iSorted],allSurvival[2,iSorted] + allErrUpper[2,iSorted],color = faceColors[2],alpha = 0.3)

    plot(tPlot*xPlotScale,out[3,:],label="|10⟩",color=edgeColors[3],lineWidth = 2);
    #errorbar(testt,allSurvival[3,:],yerr=hcat(allErrLower[3,:],allErrUpper[3,:])',color = edgeColors[3],mec=edgeColors[3],mfc=faceColors[3],linestyle="none",marker="o",capsize = 3)
    plot(testt*xPlotScale,allSurvival[3,:],color = edgeColors[3],mec=edgeColors[3],mfc=faceColors[3],linestyle="none",marker="o",markerSize = 7)
    fill_between(testt[iSorted]*xPlotScale,allSurvival[3,iSorted] - allErrLower[3,iSorted],allSurvival[3,iSorted] + allErrUpper[3,iSorted],color = faceColors[3],alpha = 0.3)

    plot(tPlot*xPlotScale,out[4,:],label="|11⟩",color=edgeColors[4],lineWidth = 2);
    #errorbar(testt,allSurvival[4,:],yerr=hcat(allErrLower[4,:],allErrUpper[4,:])',color = edgeColors[4],mec=edgeColors[4],mfc=faceColors[4],linestyle="none",marker="o",capsize = 3)
    plot(testt*xPlotScale,allSurvival[4,:],color = edgeColors[4],mec=edgeColors[4],mfc=faceColors[4],linestyle="none",marker="o",markerSize = 7)
    fill_between(testt[iSorted]*xPlotScale,allSurvival[4,iSorted] - allErrLower[4,iSorted],allSurvival[4,iSorted] + allErrUpper[4,iSorted],color = faceColors[4],alpha = 0.3)

    legend(fontsize=14);
    xlabel("Interaction time (ms)", fontsize=16);
    ylabel("Population", fontsize=16);
    ylim([0,1.01]);
    =#

end

print(solBoots)
solDf = DataFrame(solBoots',:auto)
CSV.write(dataPath*"/bootstrapParams"*dStamp*"_"*tStamp*".csv",solDf)


figure(101)
subplot(5,1,1)
hist(solDf.x1)
xlabel("Interaction time")
subplot(5,1,2)
hist(solDf.x2)
xlabel("Detuning")
subplot(5,1,3)
hist(solDf.x3)
xlabel("COM Detuning noise")
subplot(5,1,4)
hist(solDf.x4)
xlabel("REL Detuning noise")
subplot(5,1,5)
hist(solDf.x5)
xlabel("J noise")