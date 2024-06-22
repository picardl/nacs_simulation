
using QuantumOptics
using PyPlot
using Optimization
using CSV
using DataFrames
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

    dets1 = [-x[2]/2];
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
    if iter%20 == 0
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

global iter = 0;

#dStamp = "20240410"#"20240402"
#tStamp = "205742"#"200923"
dStamp = "20240419" #2 um data combined
tStamp = "190010"
#dStamp = "20240503" #2.3 um data
#tStamp = "132113"
#dStamp = "20240501" #2.6 um data
#tStamp = "10534"
#dStamp = "000000" #Davids simulated data with astigmatism = 0.2 and alpha =1
#tStamp = "000000"


dataPath = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp
#dataPath = "C:/projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

xScale = 1e-3;

(testt, survival00, errLower00, errUpper00) = load_and_extract_data(dataPath*"/00_data"*dStamp*"_"*tStamp*".csv")
(_,survival01, errLower01, errUpper01) = load_and_extract_data(dataPath*"/01_data"*dStamp*"_"*tStamp*".csv")
(_,survival10, errLower10, errUpper10) = load_and_extract_data(dataPath*"/10_data"*dStamp*"_"*tStamp*".csv")
(_,survival11, errLower11, errUpper11) = load_and_extract_data(dataPath*"/11_data"*dStamp*"_"*tStamp*".csv")

testt = testt*xScale .+ 1e-6;

#=
survival00, errLower00, errUpper00 = stirapNorm(stirap00,stirapErrLower00,stirapErrUpper00,survival00,errLower00,errUpper00);
survival01, errLower01, errUpper01 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival01, errLower01, errUpper01)
survival10, errLower10, errUpper10 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival10, errLower10, errUpper10)
survival11, errLower11, errUpper11 = stirapNorm(stirap00, stirapErrLower00,stirapErrUpper00, survival11, errLower11, errUpper11)
=#

allSurvival = hcat(survival00, survival01, survival10, survival11)';
allErrLower = hcat(errLower00, errLower01, errLower10, errLower11)' .+ 5e-3;
allErrUpper = hcat(errUpper00, errUpper01, errUpper10, errUpper11)' .+ 5e-3;

filename = "/beta0.3772-dr50.-rmax400.-astig138.32-dchi0.-dist1.79e-6-Om38000.0-del0.0-gam25.0-U2.446-nlvl25-eres200.0-tol1.0e-5-iter4000-krylov1.0e-9-dt0.0001.csv";
file_path = dataPath*filename
theory = CSV.File(file_path) |> DataFrame


edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

iSorted = sortperm(testt)
xPlotScale = 1e3;

rc("font",family="sans",size = 7)
figure(3);
plot(xPlotScale*theory.t,theory.P00,label="|00⟩",color=edgeColors[1]);

plot(testt*xPlotScale,allSurvival[1,:],color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o")
fill_between(testt[iSorted]*xPlotScale,allSurvival[1,iSorted] - allErrLower[1,iSorted],allSurvival[1,iSorted] + allErrUpper[1,iSorted],color = faceColors[1],alpha = 0.3)

plot(xPlotScale*theory.t,theory.P01,label="|0e⟩",color=edgeColors[2]);

plot(testt*xPlotScale,allSurvival[2,:],color = edgeColors[2],mec=edgeColors[2],mfc=faceColors[2],linestyle="none",marker="o")
fill_between(testt[iSorted]*xPlotScale,allSurvival[2,iSorted] - allErrLower[2,iSorted],allSurvival[2,iSorted] + allErrUpper[2,iSorted],color = faceColors[2],alpha = 0.3)

plot(xPlotScale*theory.t,theory.P10,label="|e0⟩",color=edgeColors[3]);

plot(testt*xPlotScale,allSurvival[3,:],color = edgeColors[3],mec=edgeColors[3],mfc=faceColors[3],linestyle="none",marker="o")
fill_between(testt[iSorted]*xPlotScale,allSurvival[3,iSorted] - allErrLower[3,iSorted],allSurvival[3,iSorted] + allErrUpper[3,iSorted],color = faceColors[3],alpha = 0.3)

plot(xPlotScale*theory.t,theory.P11,label="|ee⟩",color=edgeColors[4]);

plot(testt*xPlotScale,allSurvival[4,:],color = edgeColors[4],mec=edgeColors[4],mfc=faceColors[4],linestyle="none",marker="o")
fill_between(testt[iSorted]*xPlotScale,allSurvival[4,iSorted] - allErrLower[4,iSorted],allSurvival[4,iSorted] + allErrUpper[4,iSorted],color = faceColors[4],alpha = 0.3)

xlim([0,maximum(xPlotScale*testt)])

legend(fontsize=7);
xlabel("Interaction time (ms)");
ylabel("Population");
ylim([0,1.01]);
xlim([0,tPlot[end]*xPlotScale])
