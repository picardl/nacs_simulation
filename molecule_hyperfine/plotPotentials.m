const = constants();

% Coupled-channel calculation in ground state of NaCs. Default configured to
% find feshbach state and lowest 2 trap states.
    Erange = [-0.03 -0.00171];
    Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

%% simulation parameters
Nx = 4e3;
rmin = 4.5;
rmax = 1.5e2;
% Erange = [-0.0088 -0.002];
% Erange = [325200 325300]*1e9*const.h/const.hartree;

bohr2angstrom = 0.529177210903;
yscale = 29.9792458*1/const.wavenum2hartree/1e3;1e-12*const.hartree/const.h;

A1S_emo = [16501.87,4.65,10509.810,5992.060,...
    4.65431678872284,0.409201927134440,-0.0316930094322587,...
    0.118889237998245,0.272082736776197,0.206826824112031,...
    -0.147778977004416,-0.147778977004416,0.00288158643376113,...
    0.839886790201783,0.358215463426036];
X1S_emo = [16317.19,3.78,10236.048,6081.142,...
    3.77999922692013,0.665828316447219,0.159655197110690,...
    0.144813242483051,0.0964993864916689,-0.0452724688520158,...
    0.320383473138789,0.0148834980740537,-0.0180900998067171,...
    -0.388373323943793];

uterms = {'A1S','b3P','c3S','b1P','X1S','a3s'};
labels = {'A^1\Sigma','b^3\Pi','c^3\Sigma','b^1\Pi','X^1\Sigma','a^3\Sigma'};

data = load(['lib/rosario_potentials/A1S_FINAL_D']);
R = round(data(:,1),5);
R = R(:)';

allPot = zeros(length(uterms),numel(R));
allPot(1,:) = data(:,2) + Eoffs;

for i = 2:4
    data = load(['lib/rosario_potentials/' uterms{i} '_FINAL_D']);
    allPot(i,:) = data(:,2) + Eoffs;
end

allPot(5,:) = NaCsXPES(R) - NaCsXPES(1e3);
allPot(6,:) = NaCsaPES(R) - NaCsXPES(1e3);

A1S_zaharova = const.wavenum2hartree*emo(bohr2angstrom*R,A1S_emo(1),A1S_emo(2),A1S_emo(3),A1S_emo(4),A1S_emo(5),A1S_emo(6:end));
b3p_zaharova = const.wavenum2hartree*emo(bohr2angstrom*R,X1S_emo(1),X1S_emo(2),X1S_emo(3),X1S_emo(4),X1S_emo(5),X1S_emo(6:end));
bEmp = NaCsbbPES(R)  - NaCsXPES(1e3) ;
% thresh = const.wavenum2hartree*emo(1e3,X1S(1),X1S(2),X1S(3),X1S(4),X1S(5),X1S(6:end));
% 
% b3p_zaharova = b3p_zaharova - 184.68*const.wavenum2hartree - 4954.24*const.wavenum2hartree;% - abs(min(allPot(5,:)));%-const.wavenum2hartree*A1S_emo(1);%+allPot(2,end);
b3p_zaharova = NaCsb0PES(bohr2angstrom*R);
A1S_zaharova = NaCsA1PES(bohr2angstrom*R);

figure(1)
clf
plot(bohr2angstrom*R,allPot([5,6,1,2,3,4],:)*yscale,'linewidth',2)
% plot(R,allPot(1:4,:)*yscale)
% plot(R,allPot(5:6,:)*yscale)
hold on
plot(bohr2angstrom*R,A1S_zaharova*yscale,'linewidth',2,'linestyle','--','color','b')
plot(bohr2angstrom*R,b3p_zaharova*yscale,'linewidth',2,'linestyle','--','color','r')
plot(bohr2angstrom*R,bEmp*yscale,'linewidth',2,'linestyle','-.','color','r')
xlim([2,18])
% ylim([-155,400])
% ylim([-0.6e4,1.6e4])
ylim([-160,500])
legend([labels([5,6,1,2,3,4]),'A1S - Zaharova','b3P - Zaharova','nacsbpes'])
ylabel('Energy from atomic threshold [THz]')
xlabel('R [Å]')
% text(13,12500,'Na(3s) + Cs(6p)')
% text(13.5,1000,'Na(3s) + Cs(6s)')
text(13,380,'Na(3s) + Cs(6p)')
text(13.5,30,'Na(3s) + Cs(6s)')

function U = emo(R,Tdis,rref,Te,De,re,a)
    R = R;
%     U = (Tdis - De) + De.*(1 - exp(-alph(R,a,rref).*(R-re))).^2;
        U = Te + De.*(1 - exp(-alph(R,a,rref).*(R-re))).^2;

    U = U;
end

function out = alph(R,a,rref)
    out = 0;
    for i = 1:length(a)
        out = out +  a(i).*((R.^3 - rref.^3)./(R.^3 + rref.^3)).^(i-1);
    end
end
