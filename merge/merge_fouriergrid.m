clear;

%%
h = 6.626e-34;
hbar = h/(2*pi);
m = 1.66e-27*156;
omega = 2*pi*30e3;

waist1 = 0.623e-6;
waist2 = 1.064e-6;
xmin = -1e-6;
xmax = 3e-6;
Nx = 400;
x = linspace(xmin,xmax,Nx)';

V1 = h*10e6;
V2 = h*12e6;
x0 = 0e-6;
x1 = 2e-6;
t0 = 0.0;
t1 = 0.05e-3;
t2 = 0.1e-3;
Nt = 400;
t = linspace(0,2*t2,Nt);
dt = t(2)-t(1);

%%
gaussian =@(x0,waist,x) exp(-(x-x0).^2/waist.^2);

V0 = V1.*(-gaussian( x0,waist1,x ));
Vf = - V2 * gaussian(x1,waist2,x);
V = V1 * linramp(1,0,t1,t2,reshape(t,1,1,Nt)).*(-gaussian( minjerk(x0,x1,t0,t1,reshape(t,1,1,Nt)),waist1,x )) - V2 * gaussian(x1,waist2,x);
Vmat = cell2mat(reshape(arrayfun(@(i) diag(V(:,:,i)),1:Nt,'un',0),1,1,Nt));

d2_dx2 = fourier_deriv_matrix(x,2); % fourier-based derivative matrix
T = -(hbar^2/(2*m))*d2_dx2; % kinetic energy operator

H0 = T + diag(V0);
H0 = (H0 + H0')/2;

Hf = T + diag(Vf);
Hf = (Hf + Hf')/2;

H = T + Vmat;
H = H/hbar;

%%
[vecs_init,~] = eigenshuffle(H0/hbar);
[vecs_fin,~] = eigenshuffle(Hf/hbar);
if contains(version,'R2020b')
    H = pagemtimes(pagemtimes(vecs_init',H),vecs_init);
else
    H = mtimesx(mtimesx(vecs_init',H),vecs_init); % go to basis of initial eigenstates
end
H = (H + conj(permute(H,[2 1 3])))/2; % enforce hermiticity

%%
psi = zeros(Nx,Nt);
psi(end,1) = 1;
for i = 1:Nt-1
    psi(:,i+1) = expm(-1i*H(:,:,i)*dt)*psi(:,i);
    psi(:,i+1) = psi(:,i+1)./sqrt(sum(abs(psi(:,i+1)).^2)); % enforce normalization
end

% figure(1);
% for i = 1:Nt
%     clf;
%     hold on;
%     plot(x*1e6,abs(vecs_init(:,:,1)*psi(:,i)).^2)
%     plot(x*1e6,V(:,1,i)/abs(V1))
%     hold off;
%     xlabel('x (micron)')
%     ylabel('|\psi|^2');
%     drawnow()
% end
% 
% figure(2);
% clf;
% hold on;
% bar(0:(Nx-1),flipud(abs(psi(:,1)).^2),'facealpha',0.5);
% bar(0:(Nx-1),flipud(abs(vecs_fin'*vecs_init*psi(:,end)).^2),'facealpha',0.5);
% hold off;

function out = minjerk(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0)*(10*trel(ind).^3 - 15*trel(ind).^4 + 6*trel(ind).^5);
end

function out = linramp(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0).*trel(ind);
end

