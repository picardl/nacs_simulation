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
save_basis = 'aFC';
Nmax = 4;

E_td = [0,1];
allX = cell(size(E_td));

B = 800e-4;0.01*1e-4;%linspace(0,1000,1e3)*1e-4;

for i = 1:length(E_td)
    allX{i} = X1Sigma(B,save_basis,2,Nmax,mtot,const,E_td(i)*const.h*1e6,deg2rad(35));
end

allE = zeros(length(allX{1}.E),length(E_td));
for i = 1:length(E_td)
    allE(:,i) = allX{i}.E;
end

%%
figure(349)
mNa = diag(allX{1}.psi'*allX{1}.ops.i_Na_z*allX{1}.psi);
mCs = diag(allX{1}.psi'*allX{1}.ops.i_Cs_z*allX{1}.psi);
N = diag(allX{1}.psi'*allX{1}.ops.J_sq*allX{1}.psi);
mN = diag(allX{1}.psi'*allX{1}.ops.J_z*allX{1}.psi);
diagQnums = ([round(mNa*2)/2,round(mCs*2)/2,round(sqrt(N)),round(mN)]);

indInit = find(sum(diagQnums == [1.5,2.5,0,0],2) == 4);
indFinal = find(sum(diagQnums(:,1:2) == [1.5,2.5],2) == 2);
clf
for i = 1:length(allX{1}.E)
    if ismember(i,indFinal)
        plot(E_td,(allE(i,:) - allE(indInit,:))/const.h/1e6,'color','red');
    else
        plot(E_td,(allE(i,:) - allE(indInit,:))/const.h/1e6,'color',[0.1,0.1,0.1]);
    end
    hold on
end
% ylim([3470,3472])
% 
%%
figure(350)

indN0 = find(sum(diagQnums(:,1:3) == [1.5,2.5,0],2) == 3);
indN1 = find(sum(diagQnums(:,1:3) == [1.5,2.5,1],2) == 3);
indN2 = find(sum(diagQnums(:,1:3) == [1.5,2.5,2],2) == 3);
indN3 = find(sum(diagQnums(:,1:3) == [1.5,2.5,3],2) == 3);
indN4 = find(sum(diagQnums(:,1:3) == [1.5,2.5,4],2) == 3);
clf
allInds = {indN0,indN1,indN2,indN3,indN4};
for j = 1:5
    subplot(5,1,j)
    for i = 1:length(allX{1}.E)
        if ismember(i,allInds{j})
            plot(E_td,(allE(i,:) - allE(i,1) -(allE(indInit,:)-allE(indInit,1)))/const.h/1e6,'color','red');
        end
        hold on
    end
end

if 0
    %%
    XTest = allX{2};
   mNa = diag(XTest.psi'*XTest.ops.i_Na_z*XTest.psi);
    mCs = diag(XTest.psi'*XTest.ops.i_Cs_z*XTest.psi);
    N = diag(XTest.psi'*XTest.ops.N_sq*XTest.psi);
%     mN = diag(XTest.psi'*(XTest.ops.J_x + XTest.ops.J_y)*XTest.psi); 
    mN = diag(XTest.psi'*XTest.ops.J_z*XTest.psi); 
    diagQnums = ([round(mNa*2)/2,round(mCs*2)/2,round(sqrt(N)),round(mN)]);
    indFinal = find(sum(diagQnums(:,1:3) == [1.5,2.5,1],2) == 3);
    nonDiag = ([round(mNa*2)/2,round(mCs*2)/2,round(sqrt(N)),mN]);
    real(nonDiag(indFinal,:))
end