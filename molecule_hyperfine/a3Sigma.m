function out = a3Sigma

const = constants();

% vibrational solver parameters
Nx = 4e4;
rmin = 4.5; % abohr
rmax = 2e2; %1e5; % abohr
% Erange = [-0.288 -0.2879];
% Erange = [-0.332 -0.287];
% Erange = [0.044 0.044+(1e12)*const.h/const.hartree]; %-0.0224 + [-1 1]*1e-4; % energy range to search, atomic units
% Erange = [const.h*100e9/const.hartree const.h*110e9/const.hartree];
% Erange = [-10e6 100e6]*const.h/const.hartree + 325130e9*const.h/const.hartree;
% Erange = [-1e-3 325130e9*const.h/const.hartree];

Eoffs = -0.331957;
Erange = [324000 326000]*1e9*const.h/const.hartree;

%% solve for vibrational energy
rtest = linspace(rmin,rmax,1e3);
V = korek_potential(rtest,'a3S')-Eoffs;

W =@(x) interp1(rtest,V,x,'linear','extrap');
[E_vib,nodes_out,psi_r,r] = ...
    cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,const.mu_nacs/const.me,1,1);

%% output
out.r = r;
out.psi = psi_r;
out.nodes = nodes_out;
out.E = E_vib;
fn = ['../data/a_' datestr(now,'YYmmDD_HHMMSS') '.mat'];
save(fn,'out')
disp(fn);

end
