clear;

const = constants();

eta = 3; % electronic state tag
Nmax = 1;
gamma_pol =0;

% simulation parameters
Nx = 500;
rmin = 4; % abohr
rmax = 15; % abohr
Erange = -0.0224 + [-1 1]*1e-4;
mtot = [2 3 4 5];
save_basis = 'UC';

B = 100e-4;0.01*1e-4;%linspace(0,1000,1e3)*1e-4;

%% build operators in uncoupled basis
bases = {'UC','IC','FC','F1C','F2C'};

basis.UC.qnums = build_basis({'eta','s_Na','s_Cs','Lambda','i_Na','i_Cs','N','S'},...
    {3,const.s_Na,const.s_Cs,0,const.i_Na,const.i_Cs,0:Nmax,0},[0 0 0 0 1 1 1 1]);
basis.UC.ops = build_operators(basis.UC.qnums);
Nstates = size(basis.UC.qnums,1);

%% successively couple angular momenta to form coupled bases
[basis.IC,basis.UC,basis.change.UC_IC] = couple_angmom(basis.UC,'i_Na','i_Cs','I');
[basis.FC,basis.IC,basis.change.IC_FC] = couple_angmom(basis.IC,'I','N','F');
[basis.F1C,basis.UC,basis.change.UC_F1C] = couple_angmom(basis.UC,'i_Na','N','F1');
[basis.F2C,basis.UC,basis.change.UC_F2C] = couple_angmom(basis.UC,'i_Cs','N','F2');
[basis.b,basis.UC,basis.change.UC_b] = couple_angmom(basis.UC,'N','S','J');

basis.change.UC_FC = basis.change.UC_IC*basis.change.IC_FC;
basis.change.F1C_FC = (basis.change.UC_F1C')*basis.change.UC_FC;
basis.change.F2C_FC = (basis.change.UC_F2C')*basis.change.UC_FC;

%% various operators
basis.UC.ops.F_z = basis.change.UC_FC * basis.FC.ops.F_z * basis.change.UC_FC';

% spin-spin scalar coupling, i_Na.i_Cs
basis.UC.ops.i_Nadoti_Cs = op_dot(basis.UC.ops,'i_Na','i_Cs');

% spin-rotation scalar coupling, i_Na.N and i_Cs.N
basis.UC.ops.i_NadotN = op_dot(basis.UC.ops,'i_Na','N');
basis.UC.ops.i_CsdotN = op_dot(basis.UC.ops,'i_Cs','N');

% electric quadrupole coupling
j1_dot_j2_cpl =@(f,j1,j2) (f.*(f+1)-j1.*(j1+1)-j2.*(j2+1))/2;
EQ =@(F,N,I) -(3*j1_dot_j2_cpl(F,N,I).^2 + 3/2*j1_dot_j2_cpl(F,N,I) - I.*(I+1).*N.*(N+1))./(2*I.*(2*I-1).*(2*N-1).*(2*N+3));
basis.F1C.ops.EQ1 = diag(EQ(basis.F1C.qnums.F1,basis.F1C.qnums.N,basis.F1C.qnums.i_Na));
basis.F2C.ops.EQ2 = diag(EQ(basis.F2C.qnums.F2,basis.F2C.qnums.N,basis.F2C.qnums.i_Cs));

% change to UC basis
basis.UC.ops.EQ1 = basis.change.UC_F1C * basis.F1C.ops.EQ1 * basis.change.UC_F1C';
basis.UC.ops.EQ2 = basis.change.UC_F1C * basis.F2C.ops.EQ2 * basis.change.UC_F1C';

%% build hamiltonian
basis.UC.ops.H0 = const.X1Sigma.Bv * basis.UC.ops.N_sq + ...
    const.X1Sigma.c4 * basis.UC.ops.i_Nadoti_Cs + ...
    const.X1Sigma.c1 * basis.UC.ops.i_NadotN + ...
    const.X1Sigma.c2 * basis.UC.ops.i_CsdotN + ...
    const.X1Sigma.eQq1 * basis.UC.ops.EQ1 + ...
    const.X1Sigma.eQq2 * basis.UC.ops.EQ2;
basis.UC.ops.Hz0 = -const.X1Sigma.gr*const.uN*( basis.UC.ops.N_z ) + ...
    -const.X1Sigma.g1*const.uN*(1-const.X1Sigma.sigma1)*( basis.UC.ops.i_Na_z ) + ...
    -const.X1Sigma.g2*const.uN*(1-const.X1Sigma.sigma2)*( basis.UC.ops.i_Cs_z );

basis.UC.ops.H = basis.UC.ops.H0 + basis.UC.ops.Hz0.*reshape(B,1,1,[]);

%% truncate basis
% basis = rmfield(basis,{'IC','FC','F1C','F2C','b'});
% [basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.F_z,mtot);

[psi,E] = eigenshuffle(basis.UC.ops.H);
E = real(E);
% initState = find((basis.UC.qnums.m_i_Na==1.5 & basis.UC.qnums.m_i_Cs==2.5 & basis.UC.qnums.N==0));

mNa = diag(psi'*basis.UC.ops.i_Na_z*psi);
mCs = diag(psi'*basis.UC.ops.i_Cs_z*psi);
N = diag(psi'*basis.UC.ops.N_sq*psi);
mN = diag(psi'*basis.UC.ops.N_z*psi);
mF = diag(psi'*basis.UC.ops.F_z*psi);
initstate = find(round(2*mNa)/2 == 1.5 & round(2*mCs)/2 == 2.5 & round(mN) == 0 & round(N) == 0);
% finalstate = find(mNa == 1.5 & mCs == 2.5 & mF == 4 & N == 2);
finalstate = find(round(2*mNa)/2 == 1.5 & round(2*mCs)/2 == 2.5 & ((round(mF) == 4 & round(N) == 2) | (round(mF) == 5 & round(N) == 2) | (round(mF) == 3 & round(N) == 2)));
Eref = E(initstate);
E = E - Eref;

% mF = basis.UC.qnums.m_i_Na + basis.UC.qnums.m_i_Cs;
% mN = basis.UC.qnums.m_N;
figure(1)
subplot(2,2,[1,3])
scatter(mF,E/const.h/1e9 + Eref,100,'Marker','_','MarkerEdgeColor','black','Linewidth',1);
hold on
scatter(mF(initstate),E(initstate)/const.h/1e9,300,'Marker','_','MarkerEdgeColor','red','Linewidth',2);
scatter(mF(finalstate),E(finalstate)/const.h/1e9,300,'Marker','_','MarkerEdgeColor','green','Linewidth',2);
ylabel("Energy [GHz]")
ylim([-0.3,3.6])
xlabel('m_F')
xlim([-6.5,6.5])
text(-4,3.3,"N = 1",'FontSize',16)
text(-4,-0.15,"N = 0",'FontSize',16);
set(gca(),'FontSize',14)

subplot(2,2,4)
scatter(mF,E/const.h/1e6,100,'Marker','_','MarkerEdgeColor','black','Linewidth',1);
hold on
scatter(mF(initstate),E(initstate)/const.h/1e6,100,'Marker','_','MarkerEdgeColor','red','Linewidth',2);
ylim([-3.5,3.5])
ylabel("Energy [MHz]")
xlabel('m_F')
xlim([-6.5,6.5])
set(gca(),'FontSize',14)

subplot(2,2,2)
scatter(mF,E/const.h/1e6,100,'Marker','_','MarkerEdgeColor','black','Linewidth',1);
hold on
scatter(mF(finalstate),E(finalstate)/const.h/1e6,100,'Marker','_','MarkerEdgeColor','green','Linewidth',2);
ylim([3.468e3,3.478e3])
ylabel("Energy [MHz]")
xlabel('m_F')
xlim([-6.5,6.5])
set(gca(),'FontSize',14)

%%
if 0
    power_up = 1; % W
    pol_up = 'sigp';

    %% laser stuff
    % electric fields
    Efield_up = 1; %sqrt(4*const.eta0*power_up/(pi*waist^2)); % V/m

    % laser spherical tensor operators
    switch pol_up
        case 'sigp'
            p_up = [1 -1i 0]/sqrt(2);
        case 'sigm'
            p_up = [1 1i 0]/sqrt(2);
        case 'pi'
            p_up = [0 0 1];
    end
    T_up = sphten(p_up);

    E0_up = mean(E(finalstate) - E(initstate));

    %% transition dipole moments, angular momentum part
    p = -1:1; % spherical index

    T1q_c_f = sphten([0 0 1]);
    T1q_X_c = sphten([1 0 0]);

    % transition dipole moments
    rot_TDM_components = operator_matrix(@transition_dipole_case_a,...
        {basis.UC.qnums,basis.UC.qnums},{'eta','J','Omega','m_J'},p,T1q_X_c);

    % off-diagonal Hamiltonian blocks
    rot_TDM_c_f = sum((-1).^reshape(p,1,1,3) .* rot_TDM_c_f_components .* reshape(T_up(p+2),1,1,3),3);
    rot_TDM_X_c = sum((-1).^reshape(p,1,1,3) .* rot_TDM_X_c_components .* reshape(T_dn(p+2),1,1,3),3);
end