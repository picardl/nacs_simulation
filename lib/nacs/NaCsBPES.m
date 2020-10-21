function V = NaCsBPES(r)

wavenum2hartree = 4.55634011e-6;
bohr2angstrom = 0.529177210903;
r = r*bohr2angstrom;

r1 = 13.18506;

R = [3.4000
    3.6000
    3.7000
    3.8000
    3.9000
    4.0000
    4.2000
    4.4000
    4.6000
    4.8000
    5.2000
    5.6000
    6.0000
    6.4000
    6.8000
    7.2000
    7.6000
    8.0000
    8.8000
    9.6000
   10.0000
   11.2000
   12.0000
   12.8000
   14.4000
   16.0000
   17.6000
   19.2000
   20.8000
   22.4000
   24.8000];

U = 1.0e+04 *[2.830337520000000
   1.975555970000000
   1.792148270000000
   1.646332090000000
   1.578764190000000
   1.558021260000000
   1.541701340000000
   1.537002310000000
   1.539498540000000
   1.546381530000000
   1.564924620000000
   1.582421580000000
   1.598106970000000
   1.612143390000000
   1.623933090000000
   1.633288270000000
   1.640474890000000
   1.645629070000000
   1.653810470000000
   1.658987970000000
   1.662458420000000
   1.664568260000000
   1.666012900000000
   1.666867710000000
   1.667777110000000
   1.668201080000000
   1.668407150000000
   1.668512590000000
   1.668569520000000
   1.668601650000000
   1.668626960000000];


U1 =@(r) interp1(R,U,r,'spline');
U2 = @(x) 16686.547 - 5.16831e7./x.^6 - 9.42143e9./x.^8 + 8.18828e11./x.^10;

V = 0*r;
V(r<=r1) = U1(r(r<=r1));
V(r>=r1) = U2(r(r>=r1));

V = (V - 4954.24)*wavenum2hartree;

end