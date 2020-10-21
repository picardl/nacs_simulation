clear;

c = constants();

%% simulation parameters
% Nx = 4000;
% rmin = 6; 4.9;
% rmax = 30;

Nx = 2000;
rmin = 5;
rmax = 100;

rtest = linspace(rmin,rmax,1e3);

%% total hamiltonian
W =@(r) NaCscPES(r);

Wtest = W(rtest);

figure(1);
clf;
plot(rtest,Wtest);
set(gca,'xscale','log')

Erange = 0.049415 + [-1 1]*1e-4; %[0.043 0.053457];

%% call the solver
[E_out,nodes_out,err_est,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,c.mu_nacs/c.me,1,1);

%% display results
diff(E_out)*c.hartree/c.h * 1e-6
errstr(E_out/c.wavenum2hartree*(29.9792458),err_est/c.wavenum2hartree*(29.9792458))

figure(2);
clf;
for i = 1:size(psi,3)
    subplot(size(psi,3),1,i)
    plot(r,psi(:,:,i),'linewidth',2);
    set(gca,'xscale','log')
    xlim([min(r) max(r)])
    ylabel('\psi(R)')
end
xlabel('R (a_0)')


% save('NaCs_scwfn_c3Sigma.mat','r','psi','E_out','nodes_out')
