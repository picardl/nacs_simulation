function plot_c_data

c_file = 'data/c_210702_111307.mat';
a_file = 'data/a_210702_105716.mat';


a_data = load(a_file);
r_a = a_data.out.r(:);
psi_a = permute(a_data.out.psi,[2 3 1]);


c_data = load(c_file);
r_c = c_data.out.r(:);
psi_c = permute(c_data.out.psi,[2 1 3]);

psi_c = interp1(r_c,psi_c,r_a);

FCF = permute(abs(trapz(r_a,psi_a.*psi_c)).^2,[2 3 1]);

figure(1);
clf;
imagesc(0:size(FCF,2)-1,0:size(FCF,1)-1,log10(FCF))
xlabel('v_{c^3\Sigma}')
ylabel('v_{a^3\Sigma}')
colorbar

figure(2);
clf;
plot(0:size(FCF,2)-1,log10(FCF(25,:)));
xlabel('v_{c^3\Sigma}')
ylabel('log_{10}(FCF)')
xlim([0 size(FCF,2)-1])

end