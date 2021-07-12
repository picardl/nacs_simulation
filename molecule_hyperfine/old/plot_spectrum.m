function plot_spectrum

c = constants();

B = 855e-4;

basis_data = load('../data/c3Sigma_basis.mat');
basis = basis_data.basis;
% ops = basis_data.ops;

% laser hamiltonian matrix elements
tdm_data = load(['../data/feshbach_c3Sigma_TDM_' num2str(B*1e4) 'G.mat']);
feshbach_TDM = tdm_data.TDM;

tdm_data = load(['../data/atom_c3Sigma_TDM_' num2str(B*1e4) 'G.mat']);
atom_TDM = tdm_data.TDM;

Be = 0.962e9*c.h;
EatomZeeman = atom_zeeman([3 3],[1 1],863e-4) - atom_zeeman([3 3],[1 1],10e-4);
f0 = 325124.967 - EatomZeeman*1e-9/c.h;


b = [
    -0.2240
    0.2516
    0.3407
    0.2279
    0.5105
    ];

% expt parameters
power_up = 1.24e-3; % W
waist = 10e-6; % m
t = 5e-6;
Efield_up = sqrt(4*c.eta0*power_up/(pi*waist^2)); % V/m
Gamma = c.h*120e6;
wef = 1e6*c.h;
% Be = 0.952 * 1e9 * c.h;
T_up = sphten([1 0 0]);
gS = 2.0023;

H_up_mol = sum((-1).^(-1:1).*feshbach_TDM.*T_up,2)*Efield_up;
H_up_atom = sum((-1).^(-1:1).*atom_TDM.*T_up,2)*Efield_up;

xmin = 325125;
xmax = 325140;
Nx = 5e2;

xfit = linspace(xmin,xmax,Nx);
[yfit,Rsc_mol,Rsc_atom] = fitfun(b,xfit);

xscale = 1;
xplotoff = 325000;

%%
figure()
clf;
hold on;
box on;
plot(xfit/xscale - xplotoff,Rsc_mol);
hold off;
ylabel('Na + Cs survival','fontsize',8)
xlabel('Pump frequency 325XXX (GHz)','fontsize',8)


    function [y,Rscatter_mol,Rscatter_atom] = fitfun(b,x)
        
        % c = constants();
        
        alpha_Na = b(1)*1e9*c.h;
        alpha_Cs = b(2)*1e9*c.h;
        rabi_scale = b(3);
        yscale = b(4);
        yoff = b(5);
        
        H = (Be)*basis.aUC.ops.Hrot ...
            + gS*basis.aUC.ops.HZ0elecspin.*B*(863/855) + alpha_Na*basis.aUC.ops.Hhf_Na ...
            + alpha_Cs*basis.aUC.ops.Hhf_Cs + wef*basis.aUC.ops.H_OmegaDoubling;
        
        % diagonalize
        [psi,E] = eig(H);
        E = real(diag(E));
        
        rabi_mol = rabi_scale * sum(psi.*H_up_mol,1)';
        rabi_atom = rabi_scale * sum(psi.*H_up_atom,1)';
        
        flaser_Hz = 1e9*(x - f0);
        
        Rscatter_mol = (Gamma/(2*c.hbar)) * sum(2*abs(rabi_mol).^2./(2*abs(rabi_mol).^2 + 4*((flaser_Hz-f0)*c.h-E).^2 + Gamma.^2),1);
        Rscatter_atom = (Gamma/(2*c.hbar)) * sum(2*abs(rabi_atom).^2./(2*abs(rabi_atom).^2 + 4*((flaser_Hz-f0)*c.h-E).^2 + Gamma.^2),1);
        y = yscale*exp(-Rscatter_mol*t) + yoff*exp(-Rscatter_atom*t);
        
    end

end