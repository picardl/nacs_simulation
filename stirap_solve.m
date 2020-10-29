clear;

const = constants();

data = load('data/transfer_Heff_aFC_sigp_sigp.mat');

%% further eliminate off-resonant states

ic = diag(data.ops.Ic)>0;

E0 = diag(data.ops.H0);

P = diag(abs(E0)/const.h < 10e6 & ic) + data.ops.IX + data.ops.If;
Q = eye(data.Nstates) - P;

EQ = mean(diag(Q*data.ops.H0*Q));

%%