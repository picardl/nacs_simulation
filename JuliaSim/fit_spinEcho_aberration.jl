using QuadGK
using PyPlot
using Optimization
using CSV
using DataFrames
using LsqFit

function proj0(
    alpha::ComplexF64,
    lambda::Float64,
    nu::Float64,
    tau::Float64,
)
    term1 = (alpha * exp(2im*nu*tau) + 2lambda * (exp(im*nu*tau)-1) - conj(alpha))
    exponent = 0.5lambda*exp(-2im*nu*tau)*(exp(im*nu*tau)-1)^2 * term1
    return 0.5 - 0.5*real(exp(exponent))
end

function proj0_avg(
    alpha_abs::Float64,
    lambda::Float64,
    nu::Float64,
    tau::Float64
)
    res = 1/(2pi) .* quadgk(0, 2pi) do phi
        return proj0(alpha_abs*exp(im*phi), lambda, nu, tau)
    end
    return res[1]
end

amu = 1.66054e-27;
m_nacs = (132.905451933+22.9897692820)*amu;
h = 2*pi*1.054571817646156e-34
nu = 2*pi*5.2e3;
zr = pi*1.19e-6^2/1.064e-6;
U = 3.5e6*h;
alpha = 0.2;

dStamp = "20240526" 
tStamp = "220507"

dataPath = "C:/projects/nacs_simulation/JuliaSim/experimentalData/"*dStamp*"_"*tStamp

xScale = 1e-6;

(testt, survival0, errLower0, errUpper0) = load_and_extract_data(dataPath*"/0_data"*dStamp*"_"*tStamp*".csv")
(_,survival1, errLower1, errUpper1) = load_and_extract_data(dataPath*"/1_data"*dStamp*"_"*tStamp*".csv")

testt = testt*xScale;
plott = range(0,maximum(testt),100)

p0 = [0.1,nu];

@. model(x, p) = proj0_avg.(0.9416,p[1],p[2],x./2)

fit = curve_fit(model, testt, survival0, p0);
fitvals = coef(fit);

hoLength = sqrt(h/(4*pi*m_nacs*fitvals[2]))
print("Displacement:")
print(fitvals[1]*hoLength*1e9*2)
print(" nm \n")

edgeColors = [[0,113/255,187/255],[49,163,84]/255,[117,107,177]/255,[220,20,20]/255,[0,109,44]/255];
faceColors = [[177,224,255]/255,[161,217,155]/255,[188,189,220]/255,[255,142,142]/255,[44,162,95]/255];

xPlotScale = 1e6;

plot(plott*xPlotScale,model(plott,fitvals))
errorbar(testt*xPlotScale,survival0,yerr = (errLower0,errUpper0),color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",markersize = 7,capsize=3)
#plot(testt*xPlotScale,survival0,color = edgeColors[1],mec=edgeColors[1],mfc=faceColors[1],linestyle="none",marker="o",markerSize = 5)
#fill_between(testt*xPlotScale,survival0 - errLower0,survival0 + errUpper0,color = faceColors[1],alpha = 0.3)
ylim([0,0.04])
xlabel("Spin echo time (us)")
ylabel("|0⟩ population")