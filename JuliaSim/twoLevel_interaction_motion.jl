using QuantumOptics
using PyPlot
include("DynamicalDecoupling.jl")
using .DynamicalDecoupling
using Interpolations
using Statistics
using LsqFit
using Profile
using Random
using Dates
using DataFrames
using CSV

function genMotionState(bases,ns)
    """
    genMotion(bases,indices)
        generates a given 6D motional state for two molecules

    Arguments:
    - ns::Vector: Vector of motional states
    Returns:
    - ns
    """
    psi = tensor((fockstate.(bases,ns))...);;
    return psi;
end

function Is(bases)
    """
    Is(bases,indices)
        generates a tensor product of identity operators for all the bases in the list

    Arguments:
    - bases::Vector: Vector of bases
    Returns:
    - ICol
    """
    if length(bases) == 0
        return 1
    end
    ICol = tensor((identityoperator.(bases))...);
    return ICol;
end

function xOp(b,zeroPoint)
    """
    xOp(b)
        Generate position operator for a Fock basis b

    Arguments:
    - b::FockBasis: Quantum Optics Toolbox Fock Basis object
    - zeroPoint::Float Zero point motion, Sqrt(hbar/(2*m*omega)) for the harmonic oscillator

    Returns:
    - xOperator
    """
    return zeroPoint*(create(b) + destroy(b));

end

h = 6.62607015e-34;
w0 = 1.19e-6;
hbar = h/(2*pi);
eps0 = 8.8541878128e-12;
mNaCs = (132.905451933+22.9897692820)*1.66053906660e-27;
U = h*5e6;#h*1.46e6;

N = 2;
tPi = 20e-6;#32.2e-6; #Microwave pi pulse time at zero detuning
Ω = 2*pi*(1/(4*tPi));

omegarad = sqrt(4*U/(mNaCs*w0^2));

omega1s = [omegarad,1.01*omegarad,0.2*omegarad];
omega2s = 1.01068978626*[omegarad,1.01*omegarad,0.2*omegarad];
zpm1s = sqrt.(hbar./(2*mNaCs.*omega1s));
zpm2s = sqrt.(hbar./(2*mNaCs.*omega2s));

allOmegas = [omega1s;omega2s];

nMax = 2;

inter_t = 1.32e-3/2;

dets1 = [0];
dets2 = [50];

XRot1, YRot1, FreeEv1, Ps1, b1 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets1);

x1 = FockBasis(nMax);
y1 = FockBasis(nMax);
z1 = FockBasis(nMax + 1);

x2 = FockBasis(nMax);
y2 = FockBasis(nMax);
z2 = FockBasis(nMax + 1);

XRot2, YRot2, FreeEv2, Ps2, b2 = DynamicalDecoupling.genNLevelOperators(N, Ω, 2*pi*dets2);

b_col = b1⊗b2⊗x1⊗y1⊗z1⊗x2⊗y2⊗z2;
bs = [b1,b2,x1,y1,z1,x2,y2,z2];

SpSm = transition(b1,1,2) ⊗ transition(b2,2,1);
P00 = Ps1[1] ⊗ Ps2[1]⊗ Is(bs[3:8]);

θ = 0/180*pi;
R = 2e-6;
d = 4.6*3.33564e-30;
J = (d/sqrt(3))^2/(4*pi*eps0*R^3)*(1-3*cos(θ)^2);
Jprime = (d/sqrt(3))^2/(4*pi*eps0);
fudge = abs(1/(J/(2*hbar))/2*pi)/inter_t;

H_int = fudge*J/hbar/2*(SpSm⊗ Is(bs[3:8]) + dagger(SpSm⊗Is(bs[3:8])));
free_col = FreeEv1⊗FreeEv2⊗Is(bs[3:8]);
P00 = Ps1[1] ⊗ Ps2[1]⊗ Is(bs[3:8]);
P01= Ps1[1] ⊗ Ps2[2]⊗ Is(bs[3:8]);
P10= Ps1[2] ⊗ Ps2[1]⊗ Is(bs[3:8]);
P11= Ps1[2] ⊗ Ps2[2]⊗ Is(bs[3:8]);

X_col = XRot1 ⊗ Is(bs[2:8]) + identityoperator(b1)⊗XRot2⊗Is(bs[3:8]);
Y_col = YRot1 ⊗ Is(bs[2:8]) + identityoperator(b1)⊗YRot2⊗Is(bs[3:8]);

# Vanilla interaction starting in 01
psi01 = nlevelstate(b1,1) ⊗ nlevelstate(b2,2)⊗genMotionState(bs[3:8],[0,0,0,0,0,0]);
tout, psi_t = timeevolution.schroedinger(range(0,0.05,1000), psi01, H_int + free_col);
exp_val = expect(P01, psi_t);
figure(1)
plot(tout,exp_val,label="|00⟩")
xlabel("Time (s)")
ylabel("Popn in |01>")
title("Interaction starting in 01")

#Interaction with motional coupling
Deltax = xOp(bs[3],zpm1s[1])⊗Is(bs[[4,5,6,7,8]]) - Is(bs[[3,4,5]])⊗xOp(bs[6],zpm2s[1])⊗Is(bs[[7,8]]);
Deltay = Is(bs[[3]])⊗xOp(bs[4],zpm1s[2])⊗Is(bs[[5,6,7,8]]) - Is(bs[[3,4,5,6]])⊗xOp(bs[7],zpm2s[2])⊗Is(bs[[8]]);
Deltaz = Is(bs[[3,4]])⊗xOp(bs[5],zpm1s[3])⊗Is(bs[[6,7,8]]) - Is(bs[[3,4,5,6,7]])⊗xOp(bs[8],zpm2s[3])

function cpl_xOP(indX,bases)
    allOps = Vector();
    
    for i = 3:8
        if indX[i]
            push!(allOps,xOp(bases[i],1))
        else
            push!(allOps,Is(bs[i]))
        end
    end
    return allOps[1]⊗allOps[2]⊗allOps[3]⊗allOps[4]⊗allOps[5]⊗allOps[6]
end


motionalOperators = Vector();
for i = 1:5
    push!(motionalOperators,Is(bs[1:i+1])⊗(allOmegas[i]*number(bs[i+2]))⊗Is(bs[i+3:8]));
end;
push!(motionalOperators,Is(bs[1:7])⊗(allOmegas[6]*number(bs[8])));
HMot = sum(motionalOperators);
#HInterMotPartial = fudge*Jprime/hbar/2*SpSm⊗(-2/R^3*Is(bs[3:8]));
#HInterMotPartial = fudge*Jprime/hbar/2*SpSm⊗((-2/R^3*Is(bs[3:8]) + 6/R^5*Deltaz^2) + Deltay^2*(6/R^5*Is(bs[3:8]) - 45/2/R^7*Deltaz^2) + 
#Deltax*((6/R^4*Is(bs[3:8]) - 30/4/R^6*Deltaz^2) + Deltay^2*(-30/R^6*Is(bs[3:8]) + 315/2/R^8*Deltaz^2))
#+ Deltax^2*((-12/R^5*Is(bs[3:8]) + 90*Deltaz^2/R^7) + Deltay^2*(90/R^7*Is(bs[3:8]) - 630/R^9*Deltaz^2)));
HInterMotPartial = fudge*Jprime/hbar/2*SpSm⊗((-2/R^3*Is(bs[3:8]) + 6/R^5*Deltaz^2) + Deltay^2*6/R^5*Is(bs[3:8])  + 
Deltax*((6/R^4*Is(bs[3:8]) - 30/4/R^6*Deltaz^2) + Deltay^2*-30/R^6*Is(bs[3:8]))
- Deltax^2*12/R^5*Is(bs[3:8]) + Deltax^3*20/R^6*Is(bs[3:8]));
HInterMot = HInterMotPartial + dagger(HInterMotPartial);

# Interaction starting 01 with motion
#=psi01 = nlevelstate(b1,1) ⊗ nlevelstate(b2,2)⊗genMotionState(bs[3:8],[0,0,0,0,0,0]);
tout, psi_t = timeevolution.schroedinger(range(0,0.05,1000), psi01, HInterMot + free_col + HMot);
exp_val = expect(P01, psi_t);
exp_n = expect(Is(bs[1:2])⊗number(bs[3])⊗Is(bs[4:8]), psi_t);
figure(2)
plot(tout,exp_val,label="|00⟩")
xlabel("Time (s)")
ylabel("Popn in |01>")
title("Interaction starting in 01")
figure(3)
plot(tout,exp_n,label="|00⟩")
xlabel("Time (s)")
ylabel("Motional excitation number")=#

nx1 = 0;
ny1 = 0;
nz1 = 0;
nx2 = 0;
ny2 = 0;
nz2 = 0;

for nx1 = 0:1:2
    for ny1 = 0:1:2
        for nz1 = 0:1:2
            for nx2 = 0:1:2
                for ny2 = 0:1:2
                    for nz2 = 0:1:2
                        mot_state = [nx1,ny1,nz1,nx2,ny2,nz2]
                        print("Computing state $(mot_state)")
                        #Ramsey spin echo interaction starting in 00
                        psi00 = nlevelstate(b1,1) ⊗ nlevelstate(b2,1)⊗genMotionState(bs[3:8],mot_state);
                        #tWaits = vcat(range(1e-6,5e-3,60),range(20e-3,25e-3,60))
                        tWait = 50e-3;
                        ts = tPi.*[1/2]
                        probePhases = 0
                        tWaits = tWait*[1];

                        tout, psi_t = DynamicalDecoupling.DDSeq(ts,tWaits,[0],X_col,Y_col, HInterMot + free_col + HMot,psi00,1e6);
                        plotSR = 1e4;
                        plotInterval = Int64(round(1e6/plotSR));
                        indLoop = 1:plotInterval:length(psi_t);

                        pfN = Array{Any, 1}(undef, length(indLoop));
                        ψfN = Array{Any, 1}(undef, length(indLoop));

                        for j = 1:length(indLoop)
                            toutf, psif = timeevolution.schroedinger([0,tPi/2], psi_t[indLoop[j]], HInterMot + free_col + HMot + cos(probePhases)*X_col + sin(probePhases)*Y_col)
                            ψfN[j] = psif[end]
                        end

                        figure(3)
                        plot(indLoop/1e6,real(expect(P00, ψfN)),label="|00⟩")
                        #plot(indLoop/1e6,real(expect(P01, ψfN)),label="|01⟩")
                        #plot(indLoop/1e6,real(expect(P10, ψfN)),label="|10⟩")
                        #plot(indLoop/1e6,real(expect(P11, ψfN)),label="|11⟩")
                        legend()
                        xlabel("Spin echo wait time")
                        ylabel("Popn")
                        ylim([0,1.01]);

                        figure(4)
                        exp_n = expect(Is(bs[1:7])⊗number(bs[8]), ψfN);
                        exp_n0 = expect(Is(bs[1:7])⊗projector(fockstate(bs[8],0)), ψfN);
                        exp_n1 = expect(Is(bs[1:7])⊗projector(fockstate(bs[8],1)), ψfN);
                        exp_n10 = expect(Is(bs[1:4])⊗projector(fockstate(bs[5],0))⊗Is(bs[6:8]), ψfN);
                        exp_n11 = expect(Is(bs[1:4])⊗projector(fockstate(bs[5],1))⊗Is(bs[6:8]), ψfN);
                        plot(indLoop/1e6,exp_n,label="nbar")
                        plot(indLoop/1e6,exp_n0,label="n2=0")
                        plot(indLoop/1e6,exp_n1,label="n2=1")
                        plot(indLoop/1e6,exp_n10,label="n1=0")
                        plot(indLoop/1e6,exp_n11,label="n1=1")
                        xlabel("Spin echo wait time")
                        ylabel("Motional state")
                        legend()

                        # Save the data to a CSV file
                        data_dir = "data"
                        mkpath(data_dir)  # Create the 'data' directory if it doesn't exist
                        current_datetime = Dates.format(now(), "yyyymmdd_HHMMSS")
                        motString = "nxyz1_$(nx1)_$(ny1)_$(nz1)_nxyz2_$(nx2)_$(ny2)_$(nz2)_";
                        csv_filename = joinpath(data_dir, "RamseyData_$(motString)_$(current_datetime).csv")

                        # Extract the data you want to save
                        data = hcat(indLoop/1e6, real(expect(P00, ψfN)), real(expect(P01, ψfN)),real(expect(P10, ψfN)),real(expect(P11, ψfN)))
                        data_df = DataFrame(tWaits=data[:, 1], Expect_P00=data[:, 2], Expect_P01=data[:, 3],Expect_P10=data[:, 4],Expect_P11=data[:, 5])

                        # Write the data to the CSV file
                        CSV.write(csv_filename, data_df, header=["tWaits", "Expect_P00",  "Expect_P01", "Expect_P10","Expect_P11"])
                    end
                end
            end
        end
    end
end

