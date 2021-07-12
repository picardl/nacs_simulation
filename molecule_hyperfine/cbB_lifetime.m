function cbB_lifetime

const = constants();

file = '../data/cbB_210701_235922.mat';

Eoffs = const.h*((4*351725718500813 + 2*335116048808294)/6)/const.hartree;

data = load(file);

E = data.out.E + Eoffs;
nodes = data.out.nodes;
r = data.out.r;
qnums_upper = data.out.qnums;
ops = data.out.ops;
psi = data.out.psi;
W = data.out.W;

Nchn = size(qnums_upper,1);
Nstates = size(psi,3);

%%
lower_terms = {'X1S','a3S'};
Nlower = numel(lower_terms);

upper_terms = cellstr(qnums_upper.term(:,1:3));
Nupper = numel(upper_terms);

Nr = numel(r);

w3d2 = zeros(Nupper,Nr);
for i = 1:Nupper
    for j = 1:Nlower
        fn = ['lib/rosario_potentials/DM_' lower_terms{j} '_' upper_terms{i} '_Full'];
        if exist(fn,'file')
            data = load(fn);
            lower_fun_name = ['NaCs' lower_terms{j}(1) 'PES'];
            if strcmp(upper_terms{i}(1),'b')
                upper_fun_name = ['NaCsbbPES'];
            else
                upper_fun_name = ['NaCs' upper_terms{i}(1) 'PES'];
            end
            Vlower = feval(lower_fun_name,r);
            Vupper = feval(upper_fun_name,r);
            omega_diff = const.hartree*(Vupper-Vlower)/const.hbar;
            dipole = const.e*const.abohr*interp1(data(:,1),data(:,2),r,'linear','extrap');
            w3d2(i,:) = w3d2(i,:) + omega_diff.^3 .* dipole.^2;
        end
    end
end
w3d2 = sum(trapz(r,abs(psi).^2.*w3d2,2),1);
gamma = w3d2/(3*pi*const.eps0*const.hbar*const.c^3);
gamma = gamma(:);

p = reshape(trapz(r,abs(psi).^2,2),Nchn,Nstates);
[~,ind] = max(p,[],1);

gamma_chn = {};
for i = 1:Nchn
    gamma_chn{i} = gamma(ind==i);
end

figure(1);
clf;
hold on
for i = 1:Nchn
    plot(0:numel(gamma_chn{i})-1,gamma_chn{i}/(2*pi))
end
hold off

end