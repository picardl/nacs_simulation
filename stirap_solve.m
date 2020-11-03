function stirap_solve

const = constants();

data = load('data/transfer_Heff_aFC_sigp_sigp.mat');
b = data.b;



power_up = 0.05; % W
power_dn = 30e-6; % W
waist = 13e-6; % m

%%
% electric fields
Efield_up = sqrt(4*const.eta0*power_up/(pi*waist^2)); % V/m
Efield_dn = sqrt(4*const.eta0*power_dn/(pi*waist^2)); % V/m

H0 = b.ops.H0;
Hu = b.ops.Hu*Efield_up;
Hd = b.ops.Hd*Efield_dn;

%% further eliminate off-resonant states
Nf = trace(b.ops.If);
NX = trace(b.ops.IX);
Nc = trace(b.ops.Ic);
Nstates = Nf+NX+Nc;

E0 = diag(b.ops.H0);

P = diag(abs(E0)/const.h < 10e6 & (diag(b.ops.Ic)>0)) + b.ops.IX + b.ops.If;
Q = eye(Nstates) - P;

EQ = mean(diag(Q*H0*Q));

GpQ = inv(diag(EQ-H0(Q>0)));
Gp = zeros(size(H0));
Gp(diag(Q)>0,diag(Q)>0) = GpQ;

Vmp = P*(Hu+Hd)*Q;
Vpm = Q*(Hu+Hd)*P;

Vmm = P*(Hu+Hd)*P;
Vpp = Q*(Hu+Hd)*Q;

Hm = P*H0*P + Vmp*Gp*Vpm;

b.ops.P = P;


% bb.b = b;
b = truncate_basis(b,@(ops) ops.P,1);

    function Heff()
        GpQ = inv(diag(EQ-H0(Q>0)));
        Gp = zeros(size(H0));
        Gp(diag(Q)>0,diag(Q)>0) = GpQ;
        
        Vmp = P*(Hu+Hd)*Q;
        Vpm = Q*(Hu+Hd)*P;
        
        Vmm = P*(Hu+Hd)*P;
        Vpp = Q*(Hu+Hd)*Q;
        
        Hm = P*H0*P + Vmp*Gp*Vpm;
    end

end