clear;

const = constants();

eta = 3; % electronic state tag
Nmax = 1;
gamma_pol =0;

% simulation parameters
Nx = 10000;
rmin = 4; % abohr
rmax = 15; % abohr
Erange = -0.0224 + [-1 1]*1e-4;
mtot = [-6:6];
save_basis = 'aUC';
Nmax = 2;

B = 860e-4;0.01*1e-4;%linspace(0,1000,1e3)*1e-4;
E_td = 0;
basis = X1Sigma(B,save_basis,2,Nmax,mtot,const,E_td*const.h*1e6,0);

%% truncate basis
% basis = rmfield(basis,{'IC','FC','F1C','F2C','b'});
% [basis,rc_keep,Nchn] = truncate_basis(basis,@(ops) ops.F_z,mtot);

psi = basis.psi;
E = basis.E;
initState = find((basis.qnums.m_i_Na==1.5 & basis.qnums.m_i_Cs==2.5 & basis.qnums.J==0));

mNa = diag(psi'*basis.ops.i_Na_z*psi);
mCs = diag(psi'*basis.ops.i_Cs_z*psi);
N = diag(psi'*basis.ops.J_sq*psi);
mN = diag(psi'*basis.ops.J_z*psi);
mF = diag(psi'*basis.ops.F_z*psi);
initstate = find(round(2*mNa)/2 == 1.5 & round(2*mCs)/2 == 2.5 & round(mN) == 0 & round(N) == 0);
% finalstate = find(mNa == 1.5 & mCs == 2.5 & mF == 4 & N == 2);
finalstate = find(round(2*mNa)/2 == 1.5 & round(2*mCs)/2 == 2.5 & ((round(mF) == 4 & round(N) == 2) | (round(mF) == 5 & round(N) == 2) | (round(mF) == 3 & round(N) == 2)));
Eref = E(initstate);
E = E - Eref;

figure(1)
clf
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
scatter(mF,E/const.h/1e6,200,'Marker','_','MarkerEdgeColor','black','Linewidth',1);
hold on
scatter(mF(initstate),E(initstate)/const.h/1e6,200,'Marker','_','MarkerEdgeColor','red','Linewidth',2.5);
ylim([-0.75,6])
ylabel("Energy [MHz]")
xlabel('m_F')
xlim([-6.5,6.5])
set(gca(),'FontSize',14)

subplot(2,2,2)
scatter(mF,E/const.h/1e6,200,'Marker','_','MarkerEdgeColor','black','Linewidth',1);
hold on
scatter(mF(finalstate),E(finalstate)/const.h/1e6,200,'Marker','_','MarkerEdgeColor','green','Linewidth',2.5);
ylim([3.470e3,3.478e3])
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