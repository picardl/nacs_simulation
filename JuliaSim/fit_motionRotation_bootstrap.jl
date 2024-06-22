
using PyPlot
using CSV
using DataFrames
using LsqFit
#using Zygote
using ForwardDiff
using Optimization
using OptimizationOptimJL
using Distributions
using HypothesisTests
import Plots

function resid(tData,yData,errLower,errUpper,theory)
    weights = 1.0./(errLower.^2.0 + errUpper.^2.0 .+ 1e-9);
    resids = ((yData - theory).^2.0).*weights;
    return sum(resids)
end

function chisq(tData,yData,errLower,errUpper,theory)
    resids = ((yData - theory).^2.0)./theory;
    return sum(resids)
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

function genBootData(survival0,trials)
    nPoints = length(survival0)
    allSurvival = zeros(nPoints)
    allErrLower = zeros(nPoints)
    allErrUpper = zeros(nPoints)

    samples = zeros(nPoints)

    dists = Vector{Any}(undef,nPoints)
    
    for i = 1:nPoints
            dists[i] = Distributions.Binomial(trials[i],survival0[i]);
            samples[i] = rand(dists[i]);
            conf = confint(HypothesisTests.BinomialTest(samples[i],trials[i]);level=0.682689);

            allSurvival[i] = samples[i]/trials[i];
            allErrLower[i] = allSurvival[i] - conf[1];
            allErrUpper[i] =  conf[2] - allSurvival[i];
    end
    return(allSurvival,allErrLower,allErrUpper)
end

if true
# aberration_type_arr = [106.4, 117.04, 127.68, 138.32,148.96, 159.6, 170.24, 180.88, 191.52, 202.16, 212.8]; # errors for dr40, rmax500, num_lvls5 ~ 10^-3
# beta_arr = [0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75] # T = 0.56 from fitting Gabriel's code; extra heating from other processes
# gam_deph_motion_arr = 1 ./ [5e-3:5e-3:40e-3;]

aberration_type_arr = [106.4, 117.04, 127.68, 138.32,148.96, 159.6, 170.24, 180.88, 191.52, 202.16, 212.8]; # errors for dr40, rmax500, num_lvls5 ~ 10^-3
beta_arr = [0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75] # T = 0.56 from fitting Gabriel's code; extra heating from other processes
gam_deph_motion_arr = sort(vcat(1 ./ [5e-3:5e-3:40e-3;],[60,70,80,90,110,120]))

# aberration_type_arr = [117.04, 127.68, 138.32,148.96, 159.6, 170.24, 180.88, 191.52]; # errors for dr40, rmax500, num_lvls5 ~ 10^-3
# beta_arr = [0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75] # T = 0.56 from fitting Gabriel's code; extra heating from other processes
# gam_deph_motion_arr = 1 ./ [10e-3:5e-3:35e-3;]

ab3 = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
beta3 = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
gam3 = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
resids = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
residsSpinEcho= Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
residsXY8 = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
residsDrive = Array{Float64}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))


xScale = [1e-6,1e-6,1e-3];
xPlotScale = 1e6;

dStampSpinEcho = "20240526"
tStampSpinEcho = "220507"
dataPathSpinEcho = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/"*dStampSpinEcho*"_"*tStampSpinEcho
(tSpinEcho, survivalSpinEcho, errLowerSpinEcho, errUpperSpinEcho) = load_and_extract_data(dataPathSpinEcho*"/0_data"*dStampSpinEcho*"_"*tStampSpinEcho*".csv")
tSpinEcho = tSpinEcho.*xScale[1];

 df = CSV.File(dataPathSpinEcho*"/rawTrials_data"*dStampSpinEcho*"_"*tStampSpinEcho*".csv") |> DataFrame
 trialsSpinEcho = df.trials;

dStampXY8 = "20240526"
tStampXY8  = "115643"
dataPathXY8  = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/"*dStampXY8*"_"*tStampXY8
(tXY8 , survivalXY8 , errLowerXY8 , errUpperXY8 ) = load_and_extract_data(dataPathXY8*"/0_data"*dStampXY8*"_"*tStampXY8*".csv")
tXY8 = tXY8.*xScale[2];
df = CSV.File(dataPathXY8*"/rawTrials_data"*dStampXY8*"_"*tStampXY8*".csv") |> DataFrame
trialsXY8 = df.trials;

dStampDrive = "20240530"
tStampDrive  = "202419"
dataPathDrive  = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/"*dStampDrive*"_"*tStampDrive
(tDrive, survivalDrive , errLowerDrive , errUpperDrive ) = load_and_extract_data(dataPathDrive*"/0_data"*dStampDrive*"_"*tStampDrive*".csv")
tDrive = tDrive.*xScale[3];
df = CSV.File(dataPathDrive*"/rawTrials_data"*dStampDrive*"_"*tStampDrive*".csv") |> DataFrame
trialsDrive = df.trials;

# df = CSV.File(dataPath*"/rawTrials_data"*dStamp*"_"*tStamp*".csv") |> DataFrame
# s0 = df."0";
# s1 = df."1";
# trials = df.trials;

theoryDir = "C:/nilab-projects/nacs_simulation/JuliaSim/20240618-theory-for-fits2/"

# Initialize vectors to store DataFrames and corresponding numbers
dataDrive = Array{DataFrame}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
dataXY8 = Array{DataFrame}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))
dataSpinEcho = Array{DataFrame}(undef, length(beta_arr), length(aberration_type_arr), length(gam_deph_motion_arr))

# Iterate over parameter combinations
for (i, beta) in enumerate(beta_arr)
    for (j, aberration_type) in enumerate(aberration_type_arr)
        for (k, gam_deph_motion) in enumerate(gam_deph_motion_arr)

            ab3[i, j, k] = aberration_type;
            beta3[i, j, k] = beta;
            gam3[i, j, k] = gam_deph_motion;

            # Construct filename based on current parameters
            filenameSpinEcho = "echo-depth2.448-nlvl20-phi_echo1.5707963267948966-beta$beta-dr50.-rmax400.-astig$aberration_type-dchi0.-gam_motion$gam_deph_motion.csv"
            file_path = joinpath(theoryDir, filenameSpinEcho)
            dataSpinEcho[i, j, k] = CSV.File(file_path) |> DataFrame
            residsSpinEcho[i, j, k] = resid(tSpinEcho,survivalSpinEcho,errLowerSpinEcho,errUpperSpinEcho,dataSpinEcho[i, j, k].p0)
            
            # Construct filename based on current parameters
            filenameXY8 = "xy8-depth2.448-nlvl20-phi_echo1.5707963267948966-beta$beta-dr50.-rmax400.-astig$aberration_type-dchi0.-gam_motion$gam_deph_motion.csv"
            file_path = joinpath(theoryDir, filenameXY8)
            dataXY8[i, j, k] = CSV.File(file_path) |> DataFrame
            residsXY8[i, j, k] = resid(tXY8,survivalXY8,errLowerXY8,errUpperXY8,dataXY8[i, j, k].p0)

            # Construct filename based on current parameters
            filenameDrive = "drive-depth2.448-nlvl20-phi_echo1.5707963267948966-beta$beta-dr50.-rmax400.-astig$aberration_type-dchi0.-gam_motion$gam_deph_motion.csv"
            file_path = joinpath(theoryDir, filenameDrive)
            dataDrive[i, j, k] = CSV.File(file_path) |> DataFrame
            residsDrive[i, j, k] = resid(tDrive,survivalDrive,errLowerDrive,errUpperDrive,dataDrive[i, j, k].p0)
        end
    end
end

end

# Iterate over parameter combinations
for (i, beta) in enumerate(beta_arr)
    for (j, aberration_type) in enumerate(aberration_type_arr)
        for (k, gam_deph_motion) in enumerate(gam_deph_motion_arr)
            residsSpinEcho[i, j, k] = resid(tSpinEcho,survivalSpinEcho,errLowerSpinEcho,errUpperSpinEcho,dataSpinEcho[i, j, k].p0)
            
            residsXY8[i, j, k] = resid(tXY8,survivalXY8,errLowerXY8,errUpperXY8,dataXY8[i, j, k].p0)

            residsDrive[i, j, k] = resid(tDrive,survivalDrive,errLowerDrive,errUpperDrive,dataDrive[i, j, k].p0)
        end
    end
end

resids = residsSpinEcho .+ residsXY8 .+ residsDrive;

figure(1)
clf;
plot_surface(beta3[:,:,1],ab3[:,:,1],resids[:,:,1],cmap=:turbo)
thisgam = gam_deph_motion_arr[1];
title("Gamma motion = $thisgam")
xlabel("Beta")
ylabel("Astigmatism")

figure(2)
clf;
plot_surface(ab3[5,:,:],gam3[5,:,:],resids[5,:,:],cmap=:turbo)
xlabel("Astigmatism")
ylabel("Gamma")

# figure(2)
# clf;
# for i = 1:length(gam_deph_motion_arr)
#     subplot(2,length(gam_deph_motion_arr)÷2,i)
#     plot_surface(beta3[:,:,i],ab3[:,:,i],resids[:,:,i],cmap=:turbo)
#     thisgam = gam_deph_motion_arr[i];
#     title("Gamma motion = $thisgam")
#     xlabel("Beta")
#     ylabel("Astigmatism")
# end

edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

residsLin = reshape(resids,1,:);
betaLin = reshape(beta3,1,:);
abLin = reshape(ab3,1,:);
gamLin = reshape(gam3,1,:);

# m(p,t) = sum((t[4] .- (p[1].*t[1].^2 .+ p[2].*t[1] +  p[3].*t[2].^2 .+ p[4].*t[2] +  p[5].*t[3].^2 .+ p[6].*t[3] .+ p[7])).^2);
# mSquare(p,b,a,g) = (p[1].*b.^2 .+ p[2].*b +  p[3].*a.^2 .+ p[4].*a +  p[5].*g.^2 .+ p[6].*g .+ p[7]);
m(p,t) = sum((t[4] .- (p[1].*t[1].^2 .+ p[2].*t[2].^2 .+ p[3].*t[3].^2 .+ p[4].*t[1].*t[2] .+ p[5].*t[2].*t[3] .+ p[6].*t[1].*t[3] .+ p[7].*t[1] .+ p[8].*t[2] .+ p[9].*t[3] .+ p[10])).^2 ./(p[1].*t[1].^2 .+ p[2].*t[2].^2 .+ p[3].*t[3].^2 .+ p[4].*t[1].*t[2] .+ p[5].*t[2].*t[3] .+ p[6].*t[1].*t[3] .+ p[7].*t[1] .+ p[8].*t[2] .+ p[9].*t[3] .+ p[10]));
mSquare(p,t) =(p[1].*t[1].^2 .+ p[2].*t[2].^2 .+ p[3].*t[3].^2 .+ p[4].*t[1].*t[2] .+ p[5].*t[2].*t[3] .+ p[6].*t[1].*t[3] .+ p[7].*t[1] .+ p[8].*t[2] .+ p[9].*t[3] .+ p[10]);

mFit(t,p) =(p[1].*t[1].^2 .+ p[2].*t[2].^2 .+ p[3].*t[3].^2 .+ p[4].*t[1].*t[2] .+ p[5].*t[2].*t[3] .+ p[6].*t[1].*t[3] .+ p[7].*t[1] .+ p[8].*t[2] .+ p[9].*t[3] .+ p[10]);

f = OptimizationFunction(m, Optimization.AutoForwardDiff())
p0 = [60,0.01,0,-1,0,0.01,80,-0.8,0,40];
p0 = 1.1.*[561.863,  0.0166041,  0.000796758,  -5.33492,  0.00366746,  -1.02009,  351.935,  -2.84605,  -0.237114,  167.331]
prob = OptimizationProblem(f, p0,[betaLin,abLin,gamLin,residsLin])
sol = solve(prob, LBFGS())
# xlabel("Depth (MHz)")
# ylabel("Least squared resid")
# print(fit.param)

figure(3)
fit3D = mSquare(sol,[beta3,ab3,gam3])
plot_surface(beta3[:,:,1],ab3[:,:,1],fit3D[:,:,1],cmap=:turbo)
xlabel("Beta")
ylabel("Astigmatism")

f2 = OptimizationFunction(mFit, Optimization.AutoForwardDiff())
t0 = [0.7,160,100];
prob2 = OptimizationProblem(f2, t0,sol)
sol2 = solve(prob2, LBFGS())
print(sol2)

realResids = copy(resids);

NBoots = 500;
solBoots = zeros(3,NBoots)

for bootInd = 1:NBoots
    print(".")
    (survivalSpinEchoBoot,errLowerSpinEcho,errUpperSpinEcho) =  genBootData(survivalSpinEcho,trialsSpinEcho);
    (survivalXY8Boot,errLowerXY8,errUpperXY8) =  genBootData(survivalXY8,trialsXY8);
    (survivalDriveBoot,errLowerDrive,errUpperDrive) =  genBootData(survivalDrive,trialsDrive);

    for (i, beta) in enumerate(beta_arr)
        for (j, aberration_type) in enumerate(aberration_type_arr)
            for (k, gam_deph_motion) in enumerate(gam_deph_motion_arr)
                # Construct filename based on current parameters
                residsSpinEcho[i, j, k] = resid(tSpinEcho,survivalSpinEchoBoot,errLowerSpinEcho,errUpperSpinEcho,dataSpinEcho[i, j, k].p0)
                
                # Construct filename based on current parameters
                residsXY8[i, j, k] = resid(tXY8,survivalXY8Boot,errLowerXY8,errUpperXY8,dataXY8[i, j, k].p0)

                residsDrive[i, j, k] = resid(tDrive,survivalDriveBoot,errLowerDrive,errUpperDrive,dataDrive[i, j, k].p0)

            end
        end
    end
    
    resids = residsSpinEcho .+ residsXY8 .+ residsDrive;
    residsLin = reshape(resids,1,:);

    # prob = OptimizationProblem(f, p0,[betaLin,abLin,gamLin,residsLin])
    # sol = solve(prob, LBFGS())
    # prob2 = OptimizationProblem(f2, t0,sol)
    # sol2 = solve(prob2, LBFGS())
    # solBoots[:,bootInd] = sol2;

    print(minimum(resids))
    for j = 1:1
        indmin = argmin(resids);
        solBoots[:,bootInd] = solBoots[:,bootInd] + [beta_arr[indmin[1]],aberration_type_arr[indmin[2]],gam_deph_motion_arr[indmin[3]]]
        resids[indmin] = resids[indmin] + 100;
    end
    solBoots[:,bootInd] = solBoots[:,bootInd]./1;
    print("\n")
    print(minimum(resids))
    print("\n")
end

# solDf = DataFrame(solBoots',:auto)
# CSV.write(dataPath*"/bootstrapParams"*dStamp*"_"*tStamp*".csv",solDf)

solConverged = solBoots[:,solBoots[1,:] .< 2]
solConverged = solConverged[:,solConverged[3,:] .> 0]

figure(101)
subplot(1,3,1)
hist(solConverged[1,:],50,rwidth=0.9)
xlabel("Beta")
subplot(1,3,2)
hist(solConverged[2,:],50,rwidth=0.9)
xlabel("Aberration")
subplot(1,3,3)
hist(solConverged[3,:],50,rwidth=0.9)
xlabel("Gamma")

print("\nBeta mean:")
print(mean(solConverged[1,:]))
print("\nBeta standard deviation:")
print(std(solConverged[1,:]))
print("\nAberration mean:")
print(mean(solConverged[2,:]))
print("\nAberration standard deviation:")
print(std(solConverged[2,:]))
print("\nGamma mean:")
print(mean(solConverged[3,:]))
print("\nGamma standard deviation:")
print(std(solConverged[3,:]))


betaInd = argmin(abs.(mean(solConverged[1,:]) .- beta_arr));
abInd = argmin(abs.(mean(solConverged[2,:]) .- aberration_type_arr));
gamInd = argmin(abs.(mean(solConverged[3,:]) .- gam_deph_motion_arr));
# betaInd = argmin(abs.(0.5998 .- beta_arr));
# abInd = argmin(abs.(162.69 .- aberration_type_arr));
# gamInd = argmin(abs.(143.09 .- gam_deph_motion_arr));
# betaInd = argmin(abs.(0.3772 .- beta_arr));
# abInd = argmin(abs.(138.33 .- aberration_type_arr));
# gamInd = argmin(abs.(68.27 .- gam_deph_motion_arr));

print("\n")
print(realResids[betaInd, abInd, gamInd])
rc("font",family="sans",size = 7)

figure(11)
errorbar(tSpinEcho*1e6, survivalSpinEcho, (errLowerSpinEcho, errUpperSpinEcho),color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
plot(tSpinEcho*1e6,dataSpinEcho[betaInd, abInd, gamInd].p0,color = edgeColors[1])
xlabel("Spin echo time (μs)")
ylabel("|0⟩ population")

figure(12)
errorbar(tXY8*1e6, survivalXY8, (errLowerXY8, errUpperXY8),color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
plot(tXY8*1e6,dataXY8[betaInd, abInd, gamInd].p0,color = edgeColors[1])
xlabel("XY-8 pulse separation (μs)")
ylabel("|0⟩ population")

figure(14)
errorbar(tXY8*1e6, 2*abs.(0.5.-survivalXY8), (2 .*errLowerXY8, 2 .*errUpperXY8),color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
plot(tXY8*1e6,2*abs.(0.5.-dataXY8[betaInd, abInd, gamInd].p0),color = edgeColors[1])
xlabel("XY-8 pulse separation (μs)")
ylabel("Ramsey coherence")

figure(13)
errorbar(tDrive*1e3, survivalDrive, (errLowerDrive, errUpperDrive),color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",capsize = 3)
plot(tDrive*1e3,dataDrive[betaInd, abInd, gamInd].p0,color = edgeColors[1])
xlabel("Drive time (ms)")
ylabel("|0⟩ population")