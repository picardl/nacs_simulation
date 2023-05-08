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

E_td = 0.75;
B = [1,860]*1e-4;0.01*1e-4;%linspace(0,1000,1e3)*1e-4;

allXB = cell(length(B),2);

for i = 1:length(B)
    allXB{i,1} = X1Sigma(B(i),save_basis,2,Nmax,mtot,const,0,0,deg2rad(90));
    allXB{i,2} = X1Sigma(B(i),save_basis,2,Nmax,mtot,const,E_td*const.h*1e6,0,deg2rad(90));
end

allE = zeros(length(allXB{1}.E),length(B),2);
for i = 1:length(B)
    allE(:,i,1) = allXB{i,1}.E;
    allE(:,i,2) = allXB{i,2}.E;
end

%%
E0 = zeros(length(B),3);
ES = zeros(length(B),3);
for i = 1:length(B)
    [indInit_0,indFinal_0] =  findInds(allXB{i,1},[1.5,2.5,0,0],[1.5,2.5,1]);
    [indInit_1,indFinal_1] =  findInds(allXB{i,2},[1.5,2.5,0,0],[1.5,2.5,1]);
    E0(i,:) = (allXB{i,1}.E(indFinal_0) - allXB{i,1}.E(indInit_0));
    ES(i,:) = (allXB{i,2}.E(indFinal_1) - allXB{i,2}.E(indInit_1));
end

figure(359)
clf
plot(B*1e4,(E0)/const.h/1e6,'color','red');
hold on
plot(B*1e4,(ES)/const.h/1e6,'color','green');
xlabel('B Field [G]')
ylabel('Transition freq [MHz]')
legend({'Trap depth = 0','Trap depth = 0','Trap depth = 0','Trap depth = 0.75 MHz','Trap depth =  0.75  MHz','Trap depth =  0.75  MHz'})
% ylim([3470,3472])

if 0
    
    %%
E0 = zeros(length(B),length(allXB{i,1}.E));
ES = zeros(length(B),length(allXB{i,1}.E));
for i = 1:length(B)
    E0(i,:) = (allXB{i,1}.E - min(allXB{i,1}.E));
    ES(i,:) = (allXB{i,2}.E - min(allXB{i,2}.E));
end
    figure(359)
    clf
    subplot(2,1,1)
    plot(B*1e4,(E0)/const.h/1e6,'color',[0.1,0.1,0.1]);
    ylim([3470,3474])
    subplot(2,1,2)
    plot(B*1e4,(ES)/const.h/1e6,'color','green');
    xlabel('B Field [G]')
    ylabel('Transition freq [MHz]')
    ylim([3470,3474])
end

function [inds_i,inds_f] =  findInds(X,qnums_i,qnums_f)
    mNa = diag(X.psi'*X.ops.i_Na_z*X.psi);
    mCs = diag(X.psi'*X.ops.i_Cs_z*X.psi);
    N = diag(X.psi'*X.ops.J_sq*X.psi);
    mN = diag(X.psi'*X.ops.J_z*X.psi);
    diagQnums = ([round(mNa*2)/2,round(mCs*2)/2,round(sqrt(N)),round(mN)]);
    inds_i = find(sum(diagQnums == qnums_i,2) == 4);
    inds_f = find(sum(diagQnums(:,1:3) == qnums_f,2) == 3);
end

