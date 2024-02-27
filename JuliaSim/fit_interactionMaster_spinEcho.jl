
using QuantumOptics
using PyPlot
using Optimization
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling

function masterSurvivals(x,tData,params)
    h = 6.62607015e-34
    hbar = h/(2*pi)
    eps0 = 8.8541878128e-12
    N = 2
    tPi = params["tPi"];
    θ = params["θ"];
    R = params["R"];
    d = params["d"];

    Ω = 2*pi*(1/(4*tPi))*[1];

    inter_t = x[1];

    dets1 = [-x[2]/2];
    dets2 = [-x[2]/2];
    gammas1 = x[3];
    gammas2 = x[3];


    J = (d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
    fudge = 2*abs(1/(J/(2*hbar))/2*pi)/inter_t;
    gammaJ = x[4];

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

    H_noise_1B = sum([sqrt(gamma) * Ps1[i+1] for (i,gamma) in enumerate(gammas1)])⊗sum([sqrt(gamma2) * Ps2[j+1] for (j,gamma2) in enumerate(gammas2)]);
    H_noise_2B = sqrt(gammaJ)*(transition(b1,1,2) ⊗ transition(b2,2,1) + dagger(transition(b1,1,2) ⊗ transition(b2,2,1)));
    H_noise = H_noise_1B + H_noise_2B;

    psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)

    amps = [1,0,1,0,1];
    phases::Vector{Float64} = [0,0,0,0,0];

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

function masterResid(x,tData,yData,errLower,errUpper,params)

    fitData = masterSurvivals(x,tData,params);

    weights = 1.0./(errLower.^2.0 + errUpper.^2.0);
    resids = ((fitData - yData).^2.0).*weights;
    return sum(resids)

end

params = Dict();
params["tPi"] = 20e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
params["θ"] = 0/180*pi;
params["R"] = 2e-6;
params["d"] = 4.6*3.33564e-30
params["inter_t"] = 2e-3; #Interaction pi time
Delta = 2*pi*200; #Site-by-site detuning in 2*pi*Hz

tPlot = range(1e-6,10e-3,60);

guess = [2e-3,Delta,0.01*Delta,2*pi*10]; #interaction time, detuning, detuning noise, interaction noise

testt = [1e-6,1.5e-3,3e-3,7.5e-3,9e-3]
testData = [0.06 0.015 0.045 0.05 0.02;0 0.01 0.002 0.002 0.005;0 0.01 0.002 0.002 0.005;0 0.025 0.005 0.005 0.02;]./0.06;
testErrLow = [0.01 0.005 0.01 0.01 0.01;0 0.005 0.002 0.002 0.005;0 0.005 0.002 0.002 0.005;0 0.005 0.002 0.002 0.005]./0.06;
testErrHigh = [0.01 0.005 0.01 0.01 0.01;0.002 0.005 0.002 0.002 0.005;0.002 0.005 0.002 0.002 0.005;0.002 0.005 0.002 0.002 0.005]./0.06;

resid = masterResid(guess,testt,testData,testErrLow,testErrHigh,params)
print(resid);

out = masterSurvivals(guess,tPlot,params);

figure(3);
plot(tPlot,out[1,:],label="|00⟩");
plot(tPlot,out[2,:],label="|01⟩");
plot(tPlot,out[3,:],label="|10⟩");
plot(tPlot,out[4,:],label="|11⟩");
legend();
xlabel("Spin echo wait time");
ylabel("Popn");
ylim([0,1.01]);