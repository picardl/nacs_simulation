
using PyPlot
using CSV
using DataFrames
using LsqFit
#using Zygote
#using ForwardDiff
using Optimization
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

global iter = 0;

dStamp = "20240526"
tStamp = "115643"
dataPath = "C:/nilab-projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

xScale = 1e-6;
xPlotScale = 1e6;

(testt, survival0, errLower0, errUpper0) = load_and_extract_data(dataPath*"/0_data"*dStamp*"_"*tStamp*".csv")
(_,survival1, errLower1, errUpper1) = load_and_extract_data(dataPath*"/1_data"*dStamp*"_"*tStamp*".csv")

testt = testt*xScale;

df = CSV.File(dataPath*"/rawTrials_data"*dStamp*"_"*tStamp*".csv") |> DataFrame
s0 = df."0";
s1 = df."1";
trials = df.trials;

theoryDir = "C:/nilab-projects/nacs_simulation/JuliaSim/20240618-theory-for-fits/xy8-depth/"

# Get a list of all CSV files in the directory
file_list = filter(x -> occursin(r"\.csv$", x), readdir(theoryDir))

# Initialize vectors to store DataFrames and corresponding numbers
dataframes = Vector{DataFrame}(undef, length(file_list))
depths = Vector{Float64}(undef, length(file_list))
resids = Vector{Float64}(undef, length(file_list))

# Regular expression to extract XX from the filename
pattern = r"depth(\d+\.\d+)-"

# Iterate over each file, load into DataFrame and extract XX
for (i, file) in enumerate(file_list)
    # Extract XX from the filename
    match_result  = match(pattern, file)
    if match_result  !== nothing
        xx = parse(Float64, match_result.captures[1])
        depths[i] = xx
    else
        error("Could not extract XX from filename: $file")
    end
    
    # Load the CSV file into a DataFrame
    file_path = joinpath(theoryDir, file)
    dataframes[i] = CSV.File(file_path) |> DataFrame

    resids[i] = resid(testt,survival0,errLower0,errUpper0, dataframes[i].p0)
end


edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

figure(1)
plot(depths[6:14],resids[6:14])

m(t,p) = p[3].*(t.-p[2]).^2 .+ p[1];
p0 = [10, 2.4,4000];
fit = curve_fit(m, depths[6:14], resids[6:14], p0)
plot(depths[6:14],m(depths[6:14], fit.param))
xlabel("Depth (MHz)")
ylabel("Least squared resid")
print(fit.param)
legend(["Data residuals","Polynomial fit"])

print("Optimal depth:")
print(fit.param[2])
print(" MHz\n")

NBoots = 300;
solBoots = zeros(1,NBoots)

for bootInd = 1:NBoots

    (allSurvival,allErrLower,allErrUpper) =  genBootData(survival0,trials);

    for i = 1:length(depths)
        resids[i] = resid(testt,allSurvival,allErrLower,allErrUpper, dataframes[i].p0)
    end
    fit = curve_fit(m, depths[6:14], resids[6:14], p0)
    solBoots[bootInd] = fit.param[2];

    
    # figure(2)
    # plot(depths,resids)

    # figure(3)
    # plot(depths,m(depths, fit.param))

end

solDf = DataFrame(solBoots',:auto)
CSV.write(dataPath*"/bootstrapParams"*dStamp*"_"*tStamp*".csv",solDf)

solConverged = solBoots[solBoots .> 1.5]

figure(101)
hist(solConverged',50,rwidth=0.9)
xlabel("Depth")

print("Standard deviation:")
print(std(solConverged))
print(" MHz/n")