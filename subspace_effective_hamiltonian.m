function [H0eff,Veff] = subspace_effective_hamiltonian(H0,V,P)

Nstates = size(H0,1);
NP = trace(P);
Q = eye(Nstates)-P;
NQ = trace(Q);



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


end