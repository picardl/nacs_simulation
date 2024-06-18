
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
#using Plots

function masterFidelity(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 2
    tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];

    Ω = 2*pi*(1/(4*tPi))*[1];

    inter_t = abs(x[1]);

    dets1 = [-abs(x[2])/2];
    dets2 = [abs(x[2])/2];
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
    phases::Vector{Float64} = [0,0,0,0,0];

    tsSpinEcho = tPi.*DynamicalDecoupling.tFracSpinEcho
    probePhases = 0

    rhof = Array{Any, 1}(undef, length(tData));

    psitarget = 1/sqrt(2)*(nlevelstate(b1,1) ⊗ nlevelstate(b2,1) + im*nlevelstate(b1,2) ⊗ nlevelstate(b2,2))
    rhotarget = psitarget⊗dagger(psitarget)
    for i = 1:1:length(tData)
        ts = [tPi/2,tData[i]/2,tPi,tData[i]/2,tPi/2]
        tTots,rhoTots = DynamicalDecoupling.generalTwoBodyPulseSeq_master(psi00,amps,ts,phases,X_col,Y_col,int_col,H_noise,params);
        rhof[i] = rhoTots[end];
    end

    fidels = [real(expect(rho,psitarget)) for rho=rhof]
    #S = [entropy_vn(ρ)/log(2) for ρ=ρ_red]

    return fidels, rhof
end

function masterSurvivals(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 2
    #=tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];=#
    tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];

    Ω = 2*pi*(1/(4*tPi))*[1];

    inter_t = abs(x[1]);

    dets1 = [-x[2]/2];
    dets2 = [x[2]/2];
    gammaCOM = abs(x[3]);
    gammaREL = abs(x[4]);


    J = 2*(d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = abs(x[5]);

    print(fudge)

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

params = Dict();
params["tPi"] = 13.4e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
params["θ"] = 0/180*pi;
params["R"] = 2e-6;
params["d"] = 4.6*3.33564e-30
params["inter_t"] = 2e-3; #Interaction pi time
Delta = 2*pi*180; #Site-by-site detuning in 2*pi*Hz

tPlot = range(1e-6,20e-3,1000);

#sol = [1.73e-3,Delta,2*pi*1.3,2*pi*57.2]; #interaction time, detuning, detuning noise, interaction noise
#sol = [0.0015904892838368506, 385.857590919076, 14.177601026613528, 286.8948281648608,100]; #interaction time, detuning, detuning noise, interaction noise
#sol = [1.39e-3, 50, 2*pi*2, 2*pi*5,2*pi*15]; #interaction time, detuning, COM detuning noise, relative detuning noise, interaction noise
#sol = [0.0014027886901355206, 383.8577416418294, 9.1293751084133, 16.65157632171084, 87.02792075993642]
sol = [0.0013978817736142034, 9.317074755089422, 12.676602503515985, 25.378517840797727, 107.76006188798979]
#sol = [0.001394700139470014, 9.317074755089422, 11, 22, 98] #worst case

params["tData"] = tPlot;
params["yData"] = [];
params["errLower"] = [];
params["errUpper"] = [];

out = masterSurvivals(sol,tPlot,params);

edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

figure(3);
plot(tPlot,out[1,:],label="|00⟩",color=edgeColors[1]);
#errorbar(testt,allSurvival[1,:],yerr=hcat(allErrLower[1,:],allErrUpper[1,:])',color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
plot(tPlot,out[2,:],label="|01⟩",color=edgeColors[2]);
#errorbar(testt,allSurvival[2,:],yerr=hcat(allErrLower[2,:],allErrUpper[2,:])',color = edgeColors[2],mec=edgeColors[2],mfc=faceColors[2],linestyle="none",marker="o",capsize = 3)
plot(tPlot,out[3,:],label="|10⟩",color=edgeColors[3]);
#errorbar(testt,allSurvival[3,:],yerr=hcat(allErrLower[3,:],allErrUpper[3,:])',color = edgeColors[3],mec=edgeColors[3],mfc=faceColors[3],linestyle="none",marker="o",capsize = 3)
plot(tPlot,out[4,:],label="|11⟩",color=edgeColors[4]);
#errorbar(testt,allSurvival[4,:],yerr=hcat(allErrLower[4,:],allErrUpper[4,:])',color = edgeColors[4],mec=edgeColors[4],mfc=faceColors[4],linestyle="none",marker="o",capsize = 3)

legend();
xlabel("Spin echo wait time");
ylabel("Popn");
ylim([0,1.01]);

#resid = masterResid(guess,params)
#print(resid);
tPlot = range(1e-6,1.6e-3,1000);
fidels, rho = masterFidelity(sol,tPlot,params);
vmax = maximum(fidels);
imax = argmax(fidels);
print("Peak achievable fidelity: ");
print(vmax);
print("\nAt time: ");
print(tPlot[imax]);
print("\n");


figure(4);
plot(tPlot,fidels,label="Bell state fidelity");
vlines(tPlot[imax],0,1,color="red",label="max")
legend();
xlabel("Spin echo wait time");
ylabel("Fidelity");
ylim([0,1.01]);