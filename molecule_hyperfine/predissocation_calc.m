function predissocation_calc

const = constants();

data = load('../data/cbBa_210719_221910.mat');
data = data.out;

b3P_row = strcmp(cellstr(data.qnums.term),'b3P');

data2 = load('../data/a_210719_223124.mat');
data2 = data2.out;

V = const.hbar.^2./(2*const.mu_nacs*(data.r*const.abohr).^2);
% delta = mean(diff(data2.E))*const.hartree;
delta = const.h*586e6;

rshift = linspace(-0.2,0.2,1e2);

terms = cellstr(data.qnums.term);
[~,ind] = max(trapz(data.r,data.psi.^2,2),[],1);
leg = terms(ind(:));

for i = 1:numel(data.E)
    a3S_ind = abs(data.E(i)-data2.E)==min(abs(data.E(i)-data2.E));
    psi1 = data.psi(b3P_row,:,i);
    for j = 1:numel(rshift)
        psi2 = interp1(data2.r'-rshift(j),data2.psi(:,:,a3S_ind)',data.r')';
        matrix_element = trapz(data.r,psi1.*V.*psi2);
        Gamma(i,j) = 2*pi/const.hbar * matrix_element^2/delta;
    end
end

plot(rshift,Gamma')
legend(leg)
xlabel('\DeltaR (a_0)')
ylabel('\Gamma (1/s)')
set(gca,'yscale','log')
% ylim([1e-4 1e6])

end