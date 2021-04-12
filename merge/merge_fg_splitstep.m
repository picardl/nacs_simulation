clear;

%%
h = 6.626e-34;
hbar = h/(2*pi);
m = 1.66e-27*23;

waist1 = 0.623e-6;
waist2 = 1.064e-6;
xmin = -1e-6;
xmax = 4.5e-6;
Nx = 400;
x = linspace(xmin,xmax,Nx)';

V1 = h*10e6;
V2 = h*2e6;
x0 = 0e-6;
x1 = 2e-6;
t0 = 0.0;
t1 = 0.1e-3;
t2 = 0.2e-3;
Nt = 400;
t = linspace(0,2*t2,Nt);
dt = t(2)-t(1);

%%
gaussian =@(x0,waist,x) exp(-(x-x0).^2/waist.^2);

V = - V1 * linramp(1,0,t1,t2,t).*(gaussian( minjerk(x0,x1,t0,t1,t),waist1,x ))...
    - V2 * gaussian(x1,waist2,x);
V = gpuArray(V);
V_prop = exp(-1i*V*dt/hbar);
V_prop = reshape(V_prop,[Nx,1,Nt]).*eye(Nx);

fftmat = fft(eye(Nx,'gpuArray'),[],1);
ifftmat = ifft(eye(Nx,'gpuArray'),[],1);

k = fourier_fvec(x,1);
T = -(hbar^2/(2*m))*((-1i*k).^2);

T_prop = diag(exp(-1i*T*dt/(2*hbar)));
T_prop = ifftmat*T_prop*fftmat;

T_x = ifftmat*diag(T)*fftmat;

tic;
if contains(version,'R2020b')
    prop = pagemtimes(T_prop,pagemtimes(V_prop,T_prop));
else
    prop = mtimesx(T_prop,mtimesx(V_prop,T_prop));
end
toc;

[vecs_init,vals_init] = eigenshuffle(T_x + diag(V(:,1)));
vecs_init = fliplr(vecs_init);
vals_init = flipud(vals_init);
[vecs_fin,vals_fin] = eigenshuffle(T_x + diag(V(:,end)));
vecs_fin = fliplr(vecs_fin);
vals_fin = flipud(vals_fin);

%%
tic;

psi = zeros(Nx,1,Nt);
psi(:,1) = vecs_init(:,1);

psi = gpuArray(psi);
prop = gpuArray(prop);

for i = 1:Nt-1
    psi(:,1,i+1) = prop(:,:,i)*psi(:,1,i);
end
toc;

psi = gather(psi);

psi = psi./sum(abs(psi).^2,1);

psi_finbasis = mtimesx(vecs_fin',psi);

% E_expect = mtimesx(psi,'c',mtimesx(H,psi));

% P_lost = sum(abs(psi_finbasis(vals_fin>0,1,end)).^2);

figure(1);
for i = 1:Nt
    clf;
    hold on;
    plot(x*1e6,abs(psi(:,1,i)).^2)
    plot(x*1e6,V(:,i)/abs(V1))
    hold off;
    xlabel('x (micron)')
    ylabel('|\psi|^2');
    drawnow()
end

zero_crossing = find(diff(vals_fin>0));

figure(2);
clf;
hold on; box on;
bar(0:(Nx-1),abs(psi(:,1)).^2,'facealpha',0.5);
bar(0:(Nx-1),abs(vecs_fin'*vecs_init*psi(:,end)).^2,'facealpha',0.5);
plot([1 1]*zero_crossing,[0 1],'-k');
hold off;
xlim([0 zero_crossing+5])

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

