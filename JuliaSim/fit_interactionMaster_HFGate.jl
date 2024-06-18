
using QuantumOptics
using PyPlot
using Optimization
using OptimizationBBO
using OptimizationOptimJL
using CSV
using DataFrames
using FiniteDiff
#using Zygote
#using ForwardDiff
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
#using Plots

function masterSurvivals(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 3
    #=tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];=#
    tPi = params[1];##Microwave pi pulse time at zero detuning
    θ = params[2]; #θ
    R  = params[3]; #Spacing R
    d = params[4]; #dipole moment d
    tHF = params[5]; #hyperfine pulse time
    detHF = params[6]; #Detuning of |1> state from |0>

    ampMagic = 1;
    ampHF = 8.125;

    Ω = [2*pi*(1/(4*tPi)),2*pi*(1/(4*tHF))*ampMagic/ampHF];

    inter_t = abs(x[1]);#Interaction pi time

    dets1 = [-abs(x[2])/2,-abs(x[2])/2 + detHF];
    dets2 = [abs(x[2])/2,abs(x[2])/2 + detHF];
    gammaCOM = abs(x[3]);
    gammaREL = abs(x[4]);


    J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = abs(x[5]);

    _, _, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets1);
    XRots1 =  [thisΩ * (transition(b1, 1, i+1) + dagger(transition(b1, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]
    YRots1 =  [thisΩ * (-im*transition(b1, 1, i+1) + im*dagger(transition(b1, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]

    _, _, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets2);
    XRots2 =  [thisΩ * (transition(b2, 1, i+1) + dagger(transition(b2, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]
    YRots2 =  [thisΩ * (-im*transition(b2, 1, i+1) + im*dagger(transition(b2, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]

    b_col = b1 ⊗ b2;

    FreeDet1 = detHF * Ps1[1];
    FreeDet2 = detHF * Ps2[1];

    X_col = embed(b_col,1,XRots1[1]) + embed(b_col,2,XRots2[1]);
    Y_col = embed(b_col,1,YRots1[1]) + embed(b_col,2,YRots2[1]);

    X_HFcol = embed(b_col,1,XRots1[2]) + embed(b_col,2,XRots2[2]);

    free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2);
    FreeDet_col = embed(b_col,1,FreeDet1) + embed(b_col,2,FreeDet2);

    H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    P00 = Ps1[2] ⊗ Ps2[2];
    P0e= Ps1[2] ⊗ Ps2[1];
    Pe0= Ps1[1] ⊗ Ps2[2];
    Pee= Ps1[1] ⊗ Ps2[1];
    P1e= Ps1[3] ⊗ Ps2[1];
    Pe1= Ps1[1] ⊗ Ps2[3];

    P01= Ps1[2] ⊗ Ps2[3];
    P10= Ps1[3] ⊗ Ps2[2];
    P11= Ps1[3] ⊗ Ps2[3];

    int_col = free_col + H_int;

    P0 = Ps1[2] ⊗ identityoperator(b2);

    #H_noise_1B = sum([sqrt(gamma) * Ps1[i+1] for (i,gamma) in enumerate(gammas1)])⊗sum([sqrt(gamma2) * Ps2[j+1] for (j,gamma2) in enumerate(gammas2)]);
    H_noise_DeltaCOM = sqrt(gammaCOM)*(-1*P00 - P11 + Pee);
    H_noise_DeltaREL = sqrt(gammaREL/2)*(-1*P0e + Pe0 -1*P1e + Pe1);
    H_noise_2B = sqrt(gammaJ)*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    H_noise = [H_noise_DeltaCOM, H_noise_DeltaREL,H_noise_2B];

    ψ0 = nlevelstate(b1,2) ⊗ nlevelstate(b2,3)

    tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
    probePhases = 0

    rhof = Array{Any, 1}(undef, length(tData));
    
    for j = 1:1:length(tData)
        print('.')
        ts = [tHF,tData[j]/2,tPi,tData[j]/2,tHF]
        amps = [ampHF,0,ampMagic,0,ampHF]
        bDet = [1,0,0,0,1] #Whether or not to detune microwave to resonance with hyperfine
        #ts = [tHF,(tData[j] - tHF)/4,tPi,(tData[j] + tHF)/2,tPi,(tData[j] - tHF)/4,tHF]
        #ts = [tHF,tData[j]/4,tPi,2*tHF + tData[j]/2,tPi,tData[j]/4,tHF]
        #amps = [ampHF,0,ampMagic,0,ampMagic,0,ampHF]
        #bDet = [1,0,0,0,0,0,1] #Whether or not to detune microwave to resonance with hyperfine
        #ts = [tHF,(tData[j] - tHF)/4,tPi,((tData[j] + tHF)/2 - 2*tHF)/2,2*tHF,((tData[j] + tHF)/2 - 2*tHF)/2,tPi,(tData[j] - tHF)/4,tHF]
        #amps = [ampHF,0,ampMagic,0,ampHF,0,ampMagic,0,ampHF]
        #bDet = [1,0,0,0,1,0,0,0,1] #Whether or not to detune microwave to resonance with hyperfine
        #ts = [tHF,tData[j],tHF]
        #amps = [ampHF,0,1*ampHF]
        #bDet = [1,0,1] #Whether or not to detune microwave to resonance with hyperfine
        bDetInv = abs.(bDet.-1);
        phases = [0,0,pi/2,0,0]
        tTots = []
        rhoTots = []
        rho = ψ0⊗dagger(ψ0);
        tEnd = 0;
        dt = 10e-9;
    
        for i = 1:length(ts)
            Δt = (t) -> Deltas; 
            #params["len"] = times[i]
            
            tPulse = [0,ts[i]];
    
            H = int_col + bDetInv[i]*(amps[i]*X_col*cos(phases[i]) + amps[i]*Y_col*sin(phases[i])) + bDet[i]*(amps[i]*X_HFcol+FreeDet_col);
    
            tout, rho_t = timeevolution.master(tPulse, rho,  H, H_noise);
    
            tTots = vcat(tTots,tout[end] .+ tEnd);
            rhoTots = vcat(rhoTots,rho_t[end]);
    
            tEnd = tTots[end];
            rho = rho_t[end];
        end

        rhof[j] = rhoTots[end];
    end

    fitData = zeros(4,length(rhof));
    fitData[1,:] = real(expect(P00, rhof));
    fitData[2,:] = real(expect(P01, rhof));
    fitData[3,:] = real(expect(P10, rhof));
    fitData[4,:] = real(expect(P11, rhof));
    #fitData[5,:] = real(expect(P0e, rhof));
    #fitData[6,:] = real(expect(Pe0, rhof));
    #fitData[7,:] = real(expect(P1e, rhof));
    #fitData[8,:] = real(expect(Pe1, rhof));
    #fitData[9,:] = real(expect(Pee, rhof));

    return fitData
end

function scanHFTime(x,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 3
    #=tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];=#
    tPi = params[1];##Microwave pi pulse time at zero detuning
    θ = params[2]; #θ
    R  = params[3]; #Spacing R
    d = params[4]; #dipole moment d
    #tHF = params[5]; #hyperfine pulse time
    detHF = params[6]; #Detuning of |1> state from |0>

    tHFs = x[1];
    tData = x[2];

    MESol = params[11];#Master equation decoherence parameters

    ampMagic = 1;
    ampHF = 8.125;


    rhof = Array{Any, 1}(undef, length(tData));

    for (j, tHF) in enumerate(tHFs)

        Ω = [2*pi*(1/(4*tPi)),2*pi*(1/(4*tHF))*ampMagic/ampHF];

        inter_t = abs(MESol[1]);#Interaction pi time

        dets1 = [-abs(MESol[2])/2,-abs(MESol[2])/2 + detHF];
        dets2 = [abs(MESol[2])/2,abs(MESol[2])/2 + detHF];
        gammaCOM = abs(MESol[3]);
        gammaREL = abs(MESol[4]);


        J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
        fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
        gammaJ = abs(MESol[5]);

        _, _, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets1);
        XRots1 =  [thisΩ * (transition(b1, 1, i+1) + dagger(transition(b1, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]

        _, _, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets2);
        XRots2 =  [thisΩ * (transition(b2, 1, i+1) + dagger(transition(b2, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]
        b_col = b1 ⊗ b2;

        FreeDet1 = detHF * Ps1[1];
        FreeDet2 = detHF * Ps2[1];

        X_col = embed(b_col,1,XRots1[1]) + embed(b_col,2,XRots2[1]);
        X_HFcol = embed(b_col,1,XRots1[2]) + embed(b_col,2,XRots2[2]);

        free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2);
        FreeDet_col = embed(b_col,1,FreeDet1) + embed(b_col,2,FreeDet2);

        H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
        P00 = Ps1[2] ⊗ Ps2[2];
        P0e= Ps1[2] ⊗ Ps2[1];
        Pe0= Ps1[1] ⊗ Ps2[2];
        Pee= Ps1[1] ⊗ Ps2[1];
        P1e= Ps1[3] ⊗ Ps2[1];
        Pe1= Ps1[1] ⊗ Ps2[3];

        P01= Ps1[2] ⊗ Ps2[3];
        P10= Ps1[3] ⊗ Ps2[2];
        P11= Ps1[3] ⊗ Ps2[3];

        int_col = free_col + H_int;

        P0 = Ps1[2] ⊗ identityoperator(b2);

        #H_noise_1B = sum([sqrt(gamma) * Ps1[i+1] for (i,gamma) in enumerate(gammas1)])⊗sum([sqrt(gamma2) * Ps2[j+1] for (j,gamma2) in enumerate(gammas2)]);
        H_noise_DeltaCOM = sqrt(gammaCOM)*(-1*P00 + Pee);
        H_noise_DeltaREL = sqrt(gammaREL/2)*(-1*P0e + Pe0);
        H_noise_2B = sqrt(gammaJ)*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
        H_noise = [H_noise_DeltaCOM, H_noise_DeltaREL,H_noise_2B];

        ψ0 = nlevelstate(b1,2) ⊗ nlevelstate(b2,3)

            print('.')
            ts = [tHF,tData[j]/2,tPi,tData[j]/2,tHF]
            amps = [ampHF,0,ampMagic,0,ampHF]
            bDet = [1,0,0,0,1]#Whether or not to detune microwave to resonance with hyperfine
            bDetInv = abs.(bDet.-1);
            tTots = []
            rhoTots = []
            rho = ψ0⊗dagger(ψ0);
            tEnd = 0;
            dt = 10e-9;
        
            for i = 1:length(ts)
                Δt = (t) -> Deltas; 
                #params["len"] = times[i]
                
                tPulse = [0,ts[i]];
        
                H = int_col + bDetInv[i]*(amps[i]*X_col) + bDet[i]*(amps[i]*X_HFcol+FreeDet_col);
        
                tout, rho_t = timeevolution.master(tPulse, rho,  H, H_noise);
        
                tTots = vcat(tTots,tout[end] .+ tEnd);
                rhoTots = vcat(rhoTots,rho_t[end]);
        
                tEnd = tTots[end];
                rho = rho_t[end];
            end

            rhof[j] = rhoTots[end];

end

    Ω = [2*pi*(1/(4*tPi)),2*pi*(1/(4*tData[end]))*ampMagic/ampHF];

    inter_t = abs(MESol[1]);#Interaction pi time

    dets1 = [-abs(MESol[2])/2,-abs(MESol[2])/2 + detHF];
    dets2 = [abs(MESol[2])/2,abs(MESol[2])/2 + detHF];


    J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = abs(MESol[5]);

    _, _, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets1);
    XRots1 =  [thisΩ * (transition(b1, 1, i+1) + dagger(transition(b1, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]

    _, _, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, dets2);
    XRots2 =  [thisΩ * (transition(b2, 1, i+1) + dagger(transition(b2, 1, i+1))) for (i, thisΩ) in enumerate(Ω)]
    b_col = b1 ⊗ b2;

    FreeDet1 = detHF * Ps1[1];
    FreeDet2 = detHF * Ps2[1];

    X_col = embed(b_col,1,XRots1[1]) + embed(b_col,2,XRots2[1]);
    X_HFcol = embed(b_col,1,XRots1[2]) + embed(b_col,2,XRots2[2]);

    free_col = embed(b_col,1,FreeEv1) + embed(b_col,2,FreeEv2);
    FreeDet_col = embed(b_col,1,FreeDet1) + embed(b_col,2,FreeDet2);

    H_int = fudge*J/hbar/2*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    P00 = Ps1[2] ⊗ Ps2[2];
    P0e= Ps1[2] ⊗ Ps2[1];
    Pe0= Ps1[1] ⊗ Ps2[2];
    Pee= Ps1[1] ⊗ Ps2[1];
    P1e= Ps1[3] ⊗ Ps2[1];
    Pe1= Ps1[1] ⊗ Ps2[3];

    P01= Ps1[2] ⊗ Ps2[3];
    P10= Ps1[3] ⊗ Ps2[2];
    P11= Ps1[3] ⊗ Ps2[3];


    fitData = zeros(9,length(rhof));
    fitData[1,:] = real(expect(P00, rhof));
    fitData[2,:] = real(expect(P01, rhof));
    fitData[3,:] = real(expect(P10, rhof));
    fitData[4,:] = real(expect(P11, rhof));
    fitData[5,:] = real(expect(P0e, rhof));
    fitData[6,:] = real(expect(Pe0, rhof));
    fitData[7,:] = real(expect(P1e, rhof));
    fitData[8,:] = real(expect(Pe1, rhof));
    fitData[9,:] = real(expect(Pee, rhof));

    return fitData
end


function masterResid(x,params)
    tData = params[7];
    yData = params[8];
    errLower = params[9];
    errUpper = params[10];

    fitData = masterSurvivals(x,tData,params);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((fitData - yData).^2.0).*weights;
    print(sum(resids))
    print("\n")
    return sum(resids)

end

function masterHFTime(x,params)
    fitData = scanHFTime(x,params);
    return 1 .- (fitData[2,:] .- fitData[3,:]);
end

function masterResidGlob(x)
    tData = globparams[7];
    yData = globparams[8];
    errLower = globparams[9];
    errUpper = globparams[10];

    fitData = masterSurvivals(x,tData,globparams);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((fitData - yData).^2.0).*weights;
    return reshape(resids,(1,28))

end


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


#dStamp = "20240515" 
#tStamp = "111209"
dStamp = "20240519" 
tStamp = "222044"

#dataPath = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/20240225_115611"
dataPath = "C:/projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

xScale = 1e-3;

(testt, survival00, errLower00, errUpper00) = load_and_extract_data(dataPath*"/00_data"*dStamp*"_"*tStamp*".csv")
(_,survival01, errLower01, errUpper01) = load_and_extract_data(dataPath*"/01_data"*dStamp*"_"*tStamp*".csv")
(_,survival10, errLower10, errUpper10) = load_and_extract_data(dataPath*"/10_data"*dStamp*"_"*tStamp*".csv")
(_,survival11, errLower11, errUpper11) = load_and_extract_data(dataPath*"/11_data"*dStamp*"_"*tStamp*".csv")

testt = testt*xScale;



allSurvival = hcat(survival00, survival01, survival10, survival11)'
allErrLower = hcat(errLower00, errLower01, errLower10, errLower11)'
allErrUpper = hcat(errUpper00, errUpper01, errUpper10, errUpper11)'


Delta = 2*pi*50; #Site-by-site detuning in 2*pi*Hz

params = (Vector{Any}(undef,11))
params[1] = 13.4e-6;##Microwave pi pulse time at zero detuning
params[2] = 0/180*pi; #θ
params[3] = 2e-6; #Spacing R
params[4] = 4.6*3.33564e-30 #dipole moment d
params[5] = 390e-6; #HF Gate pi time
params[6] = 0#2*pi*1.89477e6;#20;#
Delta = 2*pi*50; #Site-by-site detuning in 2*pi*Hz

tPlot = range(1e-6,1.4e-3,60);
#tPlot = sort(testt);

#guess = [0.0014, 10, 10, 10, 107.48241336443607]; #interaction time, detuning, COM detuning noise, Relative detuning noise, interaction noise
#guess = [0.0014001081036026525, 0, 0, 0, 0]; #interaction time, detuning, COM detuning noise, Relative detuning noise, interaction noise
guess =[0.0013978817736142034, 9.317074755089422, 12.676602503515985, 25.378517840797727, 107.76006188798979]
#guess = [0.0012570604157658202, -419.5123134564532, 1098.088433006506, -856.7937792626279, -5.802948896867985e-5]
sol = guess

params[7] = testt;
params[8] = allSurvival;
params[9] = allErrLower;
params[10] = allErrUpper;
params[11] = sol;

paramTuple = Tuple(x for x in params)

#prob = OptimizationProblem(masterResid, guess,paramTuple, lb = [1e-3,0,0,0,0], ub = [2e-3,2*pi*100,2*pi*100,2*pi*100,2*pi*100])
prob = OptimizationProblem(masterResid, guess,paramTuple)

#prob = OptimizationProblem(masterHFTime, [300e-6,350e-6],paramTuple)
#prob = OptimizationProblem(masterHFTime, [300e-6,350e-6],paramTuple, lb = [250e-6,100e-6], ub = [2e-3,2e-3])


#sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); callback = callback_function, local_reltol = 1e-3,local_abstol = 1e-3,maxiters = 10)
#sol = solve(prob, ParticleSwarm(); callback = callback_function,x_tol = 1e-5, f_tol = 1e-3,maxiters = 5000);
#sol = solve(prob, SAMIN(rt = 0.75); callback = callback_function, x_tol = 1e-6, f_tol = 1e-3)
#sol = solve(prob, NelderMead(lower = [1e-3,-2*pi*1000,0,0,0],upper=[4e-3,2*pi*1000,2*pi*100,2*pi*100,2*pi*100]); callback = callback_function, x_tol = 1e-7, f_tol = 1e-7,maxiters = 5000)

print(sol)

#solHF = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited(); callback = callback_function, local_reltol = 1e-4,local_abstol = 1e-4,maxiters = 2000)
#solHF = solve(prob, NelderMead(lower = [200e-6,300e-6],upper=[1e-3,2e-3]); callback = callback_function, x_tol = 1e-7, f_tol = 1e-7,maxiters = 5000)

#params[5] = solHF[1];
#paramTuple = Tuple(x for x in params);
#print(solHF)

#resid = masterResid(guess,params)
#print(resid);

#=
tPlot = range(1e-6,300e-6,5);
out = testHFRabi(guess,tPlot,paramTuple);
plot(tPlot*xPlotScale,out[4,:],label="|11⟩",color=edgeColors[4],lineWidth = 2);
plot(tPlot*xPlotScale,out[5,:],label="|ee⟩",color=edgeColors[1],lineWidth = 2);
legend()
=#

out = masterSurvivals(sol,tPlot,paramTuple);
#out = masterSurvivals(guess,tPlot,paramTuple);
#out = scanHFTime(guess,tPlot,paramTuple);

edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

iSorted = sortperm(testt)
xPlotScale = 1e3;

rc("font",family="sans",size = 7)
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

#plot(tPlot*xPlotScale,out[5,:],label="|0e⟩",color=edgeColors[2],lineWidth = 2,linestyle="--");

#plot(tPlot*xPlotScale,out[6,:],label="|e0⟩",color=edgeColors[3],lineWidth = 2,linestyle="--");


legend();
xlabel("Interaction time (ms)");
ylabel("Population");
ylim([0,1.01]);


#=#Error estimation
jb = FiniteDiff.finite_difference_jacobian(masterResidGlob,sol)
weights = 1.0./(allErrLower.^2.0 + allErrUpper.^2.0 .+ 1e-9);
ww = reshape(weights,(1,116));
W = diagm(vec(ww));
cv = masterResid(sol,params)'*inv(jb'*W*jb)/(length(testt) - length(sol));
solErrs = sqrt.(diag(cv));=#