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

xmin = -0.3e-6;
xmax = 0.3e-6;
Nx = 200;

zmin = -0.5e-6;
zmax = 1.5e-6;
Nz = 200;

x = linspace(xmin,xmax,Nx)';
z = linspace(zmin,zmax,Nz)';
[X,Z] = ndgrid(x,z);


V1 = h*10e6;
V2 = h*2e6;
z0 = 0;
z1 = 1e-6;
x0 = 0e-6;
x1 = 10e-6;
t0 = 0.0;
t1 = 1e-3;
t2 = 2e-3;
Nt = 3e3;
t = reshape(linspace(0,1.5*t2,Nt),1,1,Nt);
dt = t(2)-t(1);

%%
gaussbeam =@(r0,z0,wl,w0,r,z) exp(-2*(r-r0).^2/w0^2)./(1+((z-z0)*wl/(pi*w0^2)).^2);

% V = - V1 * linramp(1,0,t1,t2,t).*(gaussbeam( minjerk(x0,x1,t0,t1,t),z0,wl1,waist1,X,Z)) ...
%     - V2 * gaussbeam(x1,z1,wl2,waist2,X,Z);

vel = -minjerk_deriv(x0,x1,t0,t1,t);

V = - V1 * linramp(1,0,t1,t2,t).*(gaussbeam( x0,z0,wl1,waist1,X,Z)) ...
    - V2 * gaussbeam(x1-minjerk(x0,x1,t0,t1,t),z1,wl2,waist2,X,Z);

V_prop = exp(-1i*(V*dt-m*vel.*X-m*vel.^2.*t)/hbar);

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

T = (K_X.^2 + K_Z.^2)*hbar^2/(2*m);
T_prop = exp(-1i*T*dt/(2*hbar));

%% find ground state of initial potential


omega_x1 = sqrt(4*V1/(m*waist1^2));
omega_z1 = sqrt(2*V1/(m*(pi*waist1^2/wl1)^2));

omega_x2 = sqrt(4*V2/(m*waist2^2));
omega_z2 = sqrt(2*V2/(m*(pi*waist2^2/wl2)^2));

aho_x1 = sqrt(hbar/(m*omega_x1));
aho_z1 = sqrt(hbar/(m*omega_z1));
aho_x2 = sqrt(hbar/(m*omega_x2));
aho_z2 = sqrt(hbar/(m*omega_z2));

Xho1 = X/aho_x1;
Zho1 = (Z-z0)/aho_z1;

Xho2 = X/aho_x2;
Zho2 = (Z-z1)/aho_z2;

psi_gs(:,:,1) = exp(-Xho1.^2/2-Zho1.^2/2);
psi_gs(:,:,2) = exp(-Xho2.^2/2-Zho2.^2/2);
psi_gs = psi_gs./sqrt(sum(sum(abs(psi_gs).^2,1),2));

%%
psi = zeros(Nx,Nz,Nt);
psi(:,:,1) = psi_gs(:,:,1);

E = zeros(Nt-1,1);
for i = 1:Nt-1
%     psi(:,:,i) = psi(:,:,i)./sqrt(sum(sum(abs(psi(:,:,i).^2),1),2));
    psi_k = fft(fft(psi(:,:,i),[],1),[],2);
    psi(:,:,i+1) = ifft(ifft(T_prop.*fft(fft(V_prop(:,:,i).*...
        ifft(ifft(T_prop.*psi_k,[],1),[],2),[],1),[],2),[],1),[],2);
    T_expect = sum(sum(conj(psi(:,:,i)).*ifft(ifft(T.*psi_k,[],1),[],2),1),2);
    V_expect = sum(sum(conj(psi(:,:,i)).*V(:,:,i).*psi(:,:,i),1),2);
    E(i) = real(T_expect + V_expect);
end

% psi = psi./sum(abs(psi).^2,1);

% psi_finbasis = mtimesx(vecs_fin',psi);

% E_expect = mtimesx(psi,'c',mtimesx(H,psi));

% P_lost = sum(abs(psi_finbasis(vals_fin>0,1,end)).^2);

gs_overlap = abs(sum(sum(conj(psi_gs(:,:,2)).*psi(:,:,end),1),2)).^2;

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
out(ind) = (x1-x0)/(t1-t0)*(30*trel(ind).^2 - 60*trel(ind).^3 + 30*trel(ind).^4);
end

function out = linramp(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0).*trel(ind);
end

function out = linramp_deriv(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
ind = trel>=0 & trel<=1;
out(ind) = (x1-x0)./(t1-t0);
end

