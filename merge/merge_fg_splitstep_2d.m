clear;

%%
h = 6.626e-34;
hbar = h/(2*pi);
m = 1.66e-27*23;

wl1 = 623e-9;
wl2 = 1064e-9;
% wl2 = 623e-9;
waist1 = wl1;
waist2 = wl2;

xmin = -1e-6;
xmax = 3e-6;
Nx = 400;

zmin = -1e-6;
zmax = 2e-6;
Nz = 400;

x = linspace(xmin,xmax,Nx)';
z = linspace(zmin,zmax,Nz)';
[X,Z] = ndgrid(x,z);


V1 = h*10e6;
V2 = h*2e6;
z0 = 0;
z1 = 1e-6;
x0 = 0e-6;
x1 = 2e-6;
t0 = 0.0;
t1 = 0.1e-3;
t2 = 1e-3;
Nt = 2e3;
t = reshape(linspace(0,2*t2,Nt),1,1,Nt);
dt = t(2)-t(1);

%%
gaussbeam =@(r0,z0,wl,w0,r,z) exp(-2*(r-r0).^2/w0^2)./(1+((z-z0)*wl/(pi*w0^2)).^2);

V = - V1 * linramp(1,0,t1,t2,t).*(gaussbeam( minjerk(x0,x1,t0,t1,t),z0,wl1,waist1,X,Z)) ...
    - V2 * gaussbeam(x1,z1,wl2,waist2,X,Z);

% V = - V1 * linramp(1,0,t1,t2,t).*(gaussbeam( x0,z0,wl1,waist1,X,Z)) ...
%     - V2 * gaussbeam(x1-minjerk(x0,x1,t0,t1,t),z1,wl2,waist2,X,Z);

V_prop = exp(-1i*V*dt/hbar);

vel = minjerk_deriv(x0,x1,t0,t1,t);

% for i = 1:10:Nt
%     pcolor(Z*1e6,X*1e6,reshape(V(:,:,i),Nx,Nz)/h*1e-6)
%     shading flat
%     axis image
%     caxis([-V1/h*1e-6 0])
%     drawnow();
% end

k_x = fourier_fvec(x,1);
k_z = fourier_fvec(z,1);
[K_X,K_Z] = ndgrid(k_x,k_z);

% K_X = K_X - m*vel/hbar.*fft(fft(ones(Nx,Nz),[],1),[],2);

T = (K_X.^2 + K_Z.^2)*hbar^2/(2*m);
T_prop = exp(-1i*T*dt/(2*hbar));

%% find ground state of initial potential
iter = 1e6;
gamma = 1e6;
alpha = 0.5;
tol = 1e-6;
T2 = T(:,:,1) - 1i*hbar*gamma;
dt2 = 1e-7;
T_prop2 = exp(-1i*T2*dt2/(2*hbar));
V2 = V(:,:,[1 Nt]) - 1i*hbar*gamma;
V_prop2 = exp(-1i*V2*dt2/(hbar));

f = 0.1;
xg = [x0 x1];
zg = [z0 z1];
w0 = [waist1 waist2];
for i = 1:2
    psi_guess(:,:,i) = exp(-(X-xg(i)).^2./(f*w0(i))^2-(Z-zg(i)).^2./(sqrt(5)*f*w0(i))^2);
end
psi_guess = psi_guess./sqrt(abs(sum(sum(psi_guess.^2,2),1)));
psi_gs = psi_guess;

figure(1); clf;
for i = 1:iter
    psi_last = psi_gs;
    for j = 1:2
        psi4 = ifft(ifft(T_prop2.*fft(fft(V_prop2(:,:,j).*ifft(ifft(T_prop2.*fft(fft(psi_gs(:,:,j),[],1),[],2),[],1),[],2),[],1),[],2),[],1),[],2);
        psi_gs(:,:,j) = (1-alpha)*psi_gs(:,:,j) + alpha*psi4;
    end
    psi_gs = psi_gs./sqrt(abs(sum(sum(psi_gs.^2,2),1)));
    err = sum(sum(sum(abs(abs(psi_gs).^2-abs(psi_last).^2)).^2));
%     subplot(2,1,1);
%     pcolor(abs(psi_gs(:,:,1)).^2);
%     shading flat;
%     axis image
%     subplot(2,1,2);
%     pcolor(abs(psi_gs(:,:,2)).^2);
%     shading flat;
%     axis image
%     drawnow();
    disp(err);
    if sum(err)<tol
        break;
    end
end

%%
psi = zeros(Nx,Nz,Nt);
psi(:,:,1) = psi_gs(:,:,1);

figure(1);
clf;
for i = 1:Nt-1
    psi(:,:,i+1) = ifft(ifft(T_prop.*fft(fft(V_prop(:,:,i).*ifft(ifft(T_prop.*fft(fft(psi(:,:,i),[],1),[],2),[],1),[],2),[],1),[],2),[],1),[],2);
end

% psi = psi./sum(abs(psi).^2,1);

% psi_finbasis = mtimesx(vecs_fin',psi);

% E_expect = mtimesx(psi,'c',mtimesx(H,psi));

% P_lost = sum(abs(psi_finbasis(vals_fin>0,1,end)).^2);

figure(1);
for i = 1:10:Nt
    clf;
    hold on;
    pcolor(z*1e6,x*1e6,abs(psi(:,:,i)).^2)
    shading flat;
    axis image
    hold off;
    xlabel('z (micron)')
    xlabel('x (micron)')
    drawnow()
end

% zero_crossing = find(diff(vals_fin>0));

% figure(2);
% clf;
% hold on; box on;
% bar(0:(Nx-1),abs(psi(:,1)).^2,'facealpha',0.5);
% bar(0:(Nx-1),abs(vecs_fin'*vecs_init*psi(:,end)).^2,'facealpha',0.5);
% plot([1 1]*zero_crossing,[0 1],'-k');
% hold off;
% xlim([0 zero_crossing+5])

function out = minjerk(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0)*(10*trel(ind).^3 - 15*trel(ind).^4 + 6*trel(ind).^5);
end

function out = minjerk_deriv(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
ind = trel>=0 & trel<=1;
out(ind) = (x1-x0)*(30*trel(ind).^2 - 60*trel(ind).^3 + 30*trel(ind).^4);
end

function out = linramp(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0).*trel(ind);
end

