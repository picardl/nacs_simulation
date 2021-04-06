clear;

x = linspace(0,2*pi,1e3);

y = cos(x);

yder = DGradient(y,x);
yder2 = DGradient( DGradient(y,x),x);

D = fourier_deriv_matrix(x,1);
D2 = fourier_deriv_matrix(x,2);

figure(1);
clf;
hold on;
% plot(x,y);
% plot(x,yder);
plot(x,yder2);
% plot(x,D*y(:),'--');
plot(x,D2*y(:),'--');