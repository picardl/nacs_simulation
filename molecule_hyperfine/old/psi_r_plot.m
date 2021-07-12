

data1=load('data/feshbach_state_800G_aUC.mat');
data2=load('data/feshbach_state_855G_aUC.mat');
data3=load('data/X1Sigma_state_855G_aUC.mat');
data4=load('data/c3Sigma_state_855G_aUC.mat');

a0 = 1e-10;

figure(1); 
clf;
hold on;
box on;

plot(data1.r/a0,sum(abs(data1.psi(:,:,1)).^2,1)*a0*2000,'linewidth',3);
plot(data2.r/a0,sum(abs(data2.psi(:,:,1)).^2,1)*a0*30,'linewidth',3);
plot([0 0],[0 0])
plot(data3.r/a0,abs(data3.psi_r).^2*a0,'linewidth',3);
plot(data4.r/a0,abs(data4.psi_r).^2*a0*3,'linewidth',3);

hold off;
set(gca,'xscale','log','fontsize',20)
xlabel('internuclear separation (A)')
ylabel('|\psi|^2')
xlim([3 30])