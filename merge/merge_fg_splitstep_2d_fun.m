function [gs_overlap,E_out,params] = merge_fg_splitstep_2d_fun(align_err,waist,t_move_ramp,power)

if nargin<4
    power = [0.4*1.8e-3 6.8e-3];
end
if nargin<3
    t_move_ramp = [0.1e-3 0.2e-3];
end

if nargin<2
    waist(1) = 0.623e-6;
    waist(2) = 1.064e-6;
end

if nargin<1
    align_err = [1e-6 0e-6];
end

params.power = power;
params.t_move_ramp = t_move_ramp;
params.waist = waist;
params.align_err = align_err;

%% constants
h = 6.626e-34;
hbar = h/(2*pi);
c = 299792458;
eps0 = 8.85418782e-12;
abohr = 5.29177210903e-11;
m = 1.66054e-27*[23 133];
alpha = [
    1531 234; % Na: 623nm, 1064nm
    -393 1163; % Cs: 623nm, 1064nm
    ] * abohr^3 * (4*pi*eps0);

%% parameters

% trap depth = alpha*E0^2/2
for i = 1:2
    for j = 1:2
        depth(i,j) = alpha(i,j) * (2*power(j)/(pi*waist(j)^2)) / (2*eps0*c);
    end
end

wl = [623e-9 1064e-9];

x0 = 0e-6;
x1 = 5e-6;
x_err = align_err(1);
z0 = 0;
z1 = align_err(2);

xmin = x0-waist(2)/2;
xmax = x_err+waist(2)/2;
Nx = 100;

zmin = z0-waist(2);
zmax = z1+waist(2);
Nz = 101;

x = linspace(xmin,xmax,Nx)';
z = linspace(zmin,zmax,Nz)';
[X,Z] = ndgrid(x,z);

t0 = 0;
t1 = t_move_ramp(1);
t2 = t1 + t_move_ramp(2);

dt = 5e-7;
t = reshape(t0:dt:1.1*t2,1,1,[]);
Nt = numel(t);

%%
gaussbeam =@(r0,z0,wl,w0,r,z) exp(-2*(r-r0).^2/w0^2)./(1+((z-z0)*wl/(pi*w0^2)).^2);

vel = -minjerk_deriv(x0,x1,t0,t1,t);

V{1} = - depth(1,1) * linramp(1,0,t1,t2,t).*(gaussbeam(x0,z0,wl(1),waist(1),X,Z)) ...
    - depth(1,2) * gaussbeam(x_err+minjerk(x1,x0,t0,t1,t),z1,wl(2),waist(2),X,Z);

V{2} = - depth(2,1) * linramp(1,0,t1,t2,t).*(gaussbeam(x_err+minjerk(-x1,x0,t0,t1,t),-z1,wl(1),waist(1),X,Z)) ...
    - depth(2,2) * gaussbeam(x0,z0,wl(2),waist(2),X,Z);

V_prop{1} = exp(-1i*(V{1}*dt-m(1)*vel.*X-m(1)*vel.^2.*t)/hbar);
V_prop{2} = exp(-1i*V{2}*dt/hbar);

% figure(1);
% clf;
% for i = 1:10:Nt
%     pcolor(Z*1e6,X*1e6,V{1}(:,:,i))
%     shading flat
%     drawnow();
% end

k_x = fourier_fvec(x,1);
k_z = fourier_fvec(z,1);
[K_X,K_Z] = ndgrid(k_x,k_z);

for i = 1:2
    T{i} = (K_X.^2 + K_Z.^2)*hbar^2/(2*m(i));
    T_prop{i} = exp(-1i*T{i}*dt/(2*hbar));
end

ind(1,:) = {[1,1],[1,2]};
ind(2,:) = {[2,2],[2,2]};

wind(1,:) = [1 2];
wind(2,:) = [2 2];

xoffs = [x0 x0+x_err; x0 x0];
zoffs = [z0 z1; z0 z0];

for j = 1:2
    for k = 1:2 % initial / final
        v0 = depth(ind{j,k}(1),ind{j,k}(2));
        
        n = wind(j,k);
        
        %% approx ground state of initial and final potential
        omega_x(j,k) = sqrt(4*v0/(m(j)*waist(n)^2));
        omega_z(j,k) = sqrt(2*v0/(m(j)*(pi*waist(n)^2/wl(n))^2));
        
        aho_x(k) = sqrt(hbar/(m(j)*omega_x(j,k)));
        aho_z(k) = sqrt(hbar/(m(j)*omega_z(j,k)));
        
        Xho{k} = (X-xoffs(j,k))/aho_x(k);
        Zho{k} = (Z-zoffs(j,k))/aho_z(k);
        
        psi_gs(:,:,k) = exp(-Xho{k}.^2/2-Zho{k}.^2/2);
    end
    psi_gs = psi_gs./sqrt(sum(sum(abs(psi_gs).^2,1),2));
    
    %%
    psi = zeros(Nx,Nz,Nt);
    psi(:,:,1) = psi_gs(:,:,1);
    
    figure(1);
    clf;
    
    E = zeros(Nt-1,1);
    for i = 1:Nt-1
        psi_k = fft(fft(psi(:,:,i),[],1),[],2);
        psi(:,:,i+1) = ifft(ifft(T_prop{j}.*fft(fft(V_prop{j}(:,:,i).*...
            ifft(ifft(T_prop{j}.*psi_k,[],1),[],2),[],1),[],2),[],1),[],2);
        T_expect = sum(sum(conj(psi(:,:,i)).*ifft(ifft(T{j}.*psi_k,[],1),[],2),1),2);
        V_expect = sum(sum(conj(psi(:,:,i)).*V{j}(:,:,i).*psi(:,:,i),1),2);
        E(i) = real(T_expect + V_expect);
        
        if ~mod(i,100)
            pcolor(Z,X,abs(psi(:,:,i)).^2);
            shading flat
            drawnow();
        end
    end
    
    E = (E(:)+depth(ind{j,2}(1),ind{j,2}(2)))/hbar/(omega_x(j,2)/2 + omega_z(j,2)/2);
    E_out(j) = E(end);
    gs_overlap(j) = abs(sum(sum(conj(psi_gs(:,:,2)).*psi(:,:,end),1),2)).^2;
end

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

