const = constants();

% Coupled-channel calculation in ground state of NaCs. Default configured to
% find feshbach state and lowest 2 trap states.
    Erange = [-0.03 -0.00171];
    Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

%% simulation parameters
Nx = 4e3;
rmin = 4.5;
rmax = 1.5e2;

bohr2angstrom = 0.529177210903;
yscale = 29.9792458*1/const.wavenum2hartree/1e3;1e-12*const.hartree/const.h;

labels = {'b^3\Pi_1 / (2)\Omega=1','c^3\Sigma^+_1 / (3)\Omega=1','b^1\Pi_1 / (4)\Omega=1','X^1\Sigma^+','a^3\Sigma^+'};

data = load(['lib/rosario_potentials/A1S_FINAL_D']);
R = round(data(:,1),5);
R = R(:)';

allPot = zeros(5,numel(R));
allPot(1,:) = data(:,2) + Eoffs;

allPot(1,:) = NaCsbbPES(R) - NaCsXPES(1e3);
allPot(2,:) = NaCscPES(R) - NaCsXPES(1e3);
allPot(3,:) = NaCsBPES(R) - NaCsXPES(1e3);
allPot(4,:) = NaCsXPES(R) - NaCsXPES(1e3);
allPot(5,:) = NaCsaPES(R) - NaCsXPES(1e3);


lc = lines(7);
ord = [4,5,1,2,3];
lc = lc([5,4,3,1,2],:);

figure(1)
clf
for i = 1:5
    plot(bohr2angstrom*R,allPot(ord(i),:)*yscale,'linewidth',2,'color',lc(i,:))
    hold on
end
xlim([2,18])
ylim([-160,500])
legend(labels(ord))
ylabel('Energy from atomic threshold [THz]')
xlabel('R [Å]')
text(12.5,380,'Na(3s) + Cs(6p_{3/2})')
text(12.5,300,'Na(3s) + Cs(6p_{1/2})')
text(12.5,30,'Na(3s) + Cs(6s)')


