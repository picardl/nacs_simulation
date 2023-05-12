function V = NaCsVd(r)
%Magnetic dip-dip and second order SO coupling for NaCs
%see http://arxiv.org/abs/2203.10362

c = constants();

wavenum2hartree = 4.556335252911937e-6;
bohr2angstrom = 0.529177210903;

Ashort = -27.8;
Along = -0.027;
Bshort = 0.8;
Blong = 0.28;

V = c.hartree*(c.alpha).^2*(Ashort.*exp(-Bshort.*r) + Along.*exp(-Blong.*r) + c.gs_Cs*c.gs_Na/(4*r.^2));

end

