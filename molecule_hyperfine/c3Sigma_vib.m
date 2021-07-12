function out = c3Sigma_vib

const = constants();

% vibrational solver parameters
Nx = 2e3;
rmin = 4.5; % abohr
rmax = 1e3; % abohr
% Erange = -0.0224 + [-1 1]*1e-4; % energy range to search, atomic units
Erange = [-1e-2 -1e-4] + NaCscPES(rmax);

%% trap
waist = 1064e-9; % meter
Power = 10e-3; % watt
I0 = 2*Power/(pi*waist^2); % W/m^2
trap_depth = (const.apol_Cs*I0)/(2*const.eps0*const.c)/const.hartree;
atrap = 2*((const.abohr)/waist).^2;
Vtrap =@(r) trap_depth*(1-exp(-atrap*r.^2));

%% solve for vibrational energy
W =@(x) NaCscPES(x) + Vtrap(x);
[E_vib,nodes_out,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,const.mu_nacs/const.me,1,1);

%% output
out.r = r;
out.psi = psi_r;
out.nodes = nodes_out;
out.E = E_vib;
fn = ['data/c_vib_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
save(fn,'out')
disp(fn)

end
