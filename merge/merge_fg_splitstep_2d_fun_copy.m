function [gs_overlap,E_out,params] = merge_fg_splitstep_2d_fun_copy(align_err,waist,t_move_ramp,power)

if nargin<4
    power = [0.3*1.8e-3 6.8e-3];
end
if nargin<3
    t_move_ramp = [0.1e-3 0.1e-3];
end

if nargin<2
    waist(1) = 0.623e-6;
    waist(2) = 1.064e-6;
end

if nargin<1
    align_err = [0e-6 0e-6];
end

save_gif = false;
animate = 2;

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
x1 = 6e-6;
x_err = align_err(1);
z0 = 0;
z1 = align_err(2);

dx = waist(2)/60;
dz = waist(2)/60;
xmax = [sqrt((waist(2)/2)^2 + x_err^2) waist(2)/2];
zmax = [sqrt(waist(2)^2 + z1^2) waist(2)];

t0 = 0;
t1 = t_move_ramp(1);
t2 = t1 + t_move_ramp(2);

dt = 5e-7;
t = t0:dt:2*t2;
tt = reshape(t,1,1,[]);
Nt = numel(tt);

%% define tweezer trajectories
tweezer{1}.x =@(t) minjerk(x0,x1,t0,t1,t);
tweezer{1}.z =@(t) minjerk(z0,z0,t0,t1,t);
tweezer{1}.I =@(t) linramp(1,0,t1,t2,t);
% tweezer{1}.I =@(t) linramp(1,1,t1,t2,t);

tweezer{2}.x =@(t) (x1+x_err)*ones(size(t));
tweezer{2}.z =@(t) z1*ones(size(t));
tweezer{2}.I =@(t) ones(size(t));

atom{1}.x =@(t) minjerk(x0,x1,t0,t1,t) - x1 + minjerk(x1,x1+x_err,t1,t2,t);
atom{1}.vx =@(t) minjerk_deriv(x0,x1,t0,t1,t) + minjerk_deriv(x1,x1+x_err,t1,t2,t);
atom{1}.ax =@(t) minjerk_accel(x0,x1,t0,t1,t) + minjerk_accel(x1,x1+x_err,t1,t2,t);

atom{1}.z =@(t) minjerk(z0,z1,t1,t2,t);
atom{1}.vz =@(t) minjerk_deriv(z0,z1,t1,t2,t);
atom{1}.az =@(t) minjerk_accel(z0,z1,t1,t2,t);

atom{2}.x =@(t) (x1+x_err)*ones(size(t));
atom{2}.z =@(t) z1*ones(size(t));
atom{2}.vx =@(t) zeros(size(t));
atom{2}.ax =@(t) zeros(size(t));
atom{2}.vz =@(t) zeros(size(t));
atom{2}.az =@(t) zeros(size(t));

gaussbeam =@(r0,z0,wl,w0,r,z) exp(-2*(r-r0).^2/w0^2)./(1+((z-z0)*wl/(pi*w0^2)).^2);

% figure(2);
% clf;
% t = t(:);
% hold on
% plot(t,DGradient(DGradient(atom{1}.x(t),t),t));
% plot(t,atom{1}.ax(t),'--');
% hold off

%%
ind(1,:) = {[1,1],[1,2]};
ind(2,:) = {[2,2],[2,2]};

wind(1,:) = [1 2];
wind(2,:) = [2 2];

xoffs = [x0 x0+x_err; x0 x0];
zoffs = [z0 z1; z0 z0];

for j = 1:2
    
    x = -xmax(j):dx:xmax(j);
    z = -zmax(j):dz:zmax(j);
    Nx = numel(x);
    Nz = numel(z);
%     z = linspace(-zmax(j),zmax(j),Nz)';
    [X0,Z0] = ndgrid(x,z);
    
    k_x = fourier_fvec(x,1);
    k_z = fourier_fvec(z,1);
    [K_X0,K_Z0] = ndgrid(k_x,k_z);
    
    first_save = true;
    
    X = X0 + atom{j}.x(tt);
    Z = Z0 + atom{j}.z(tt);
    
    V = zeros(Nx,Nz,Nt);
    psi_gs = zeros(Nx,Nz,2);
    for k = 1:2 % initial / final tweezer
        v0 = depth(ind{j,k}(1),ind{j,k}(2));
%         
        n = wind(j,k);
%         
%         %% approx ground state of initial and final potential
        omega_x(j,k) = sqrt(4*v0/(m(j)*waist(n)^2));
        omega_z(j,k) = sqrt(2*v0/(m(j)*(pi*waist(n)^2/wl(n))^2));
        
        aho_x(k) = sqrt(hbar/(m(j)*omega_x(j,k)));
        aho_z(k) = sqrt(hbar/(m(j)*omega_z(j,k)));
        
        Xho = X0/aho_x(k);
        Zho = Z0/aho_z(k);
        
        psi_gs(:,:,k) = exp(-Xho.^2/2-Zho.^2/2);
        
        %% potential
        V = V - depth(j,k) * tweezer{k}.I(tt).*...
            gaussbeam(tweezer{k}.x(tt), tweezer{k}.z(tt), ...
            wl(k), waist(k), X, Z);
        
    end
    psi_gs = psi_gs./sqrt(sum(sum(abs(psi_gs).^2)));
    
    % inertial potential (i.e. fictitious force)
    V_inert = m(j)*(atom{j}.ax(tt).*X+atom{j}.az(tt).*Z); 
    V_prop = exp(-1i*dt/hbar*(V+V_inert));
    
    T = (K_X0.^2 + K_Z0.^2)*hbar^2/(2*m(j));
    T_prop = exp(-1i*T*dt/(2*hbar));
    
    %%
%     figure(1);
%     clf;
%     for k = 1:10:numel(t)
%         subplot(1,2,1);
%         contourf(Z(:,:,k),X(:,:,k),V(:,:,k)/depth(1,1))
%         caxis([-1 0])
%         colorbar
%         
%         subplot(1,2,2);
%         contourf(Z(:,:,k),X(:,:,k),V_inert(:,:,k)/depth(1,1),10);
% %         caxis([-1 1]*20)
%         colorbar
%         
%         drawnow();
%     end
    
    %%
    psi = zeros(Nx,Nz,Nt);
    psi(:,:,1) = psi_gs(:,:,1);
    
    E = zeros(Nt-1,1);
    for i = 1:Nt-1
        psi_k = fft(fft(psi(:,:,i),[],1),[],2);
%                 psi(:,:,i+1) = ifft(ifft(T_prop(:,:,i).*fft(fft(V_prop(:,:,i).*...
%                     ifft(ifft(T_prop(:,:,i).*psi_k,[],1),[],2),[],1),[],2),[],1),[],2);
%                 T_expect = sum(sum(conj(psi(:,:,i)).*ifft(ifft(T(:,:,i).*psi_k,[],1),[],2),1),2);
                psi(:,:,i+1) = ifft(ifft(T_prop.*fft(fft(V_prop(:,:,i).*...
                    ifft(ifft(T_prop.*psi_k,[],1),[],2),[],1),[],2),[],1),[],2);
                T_expect = sum(sum(conj(psi(:,:,i)).*ifft(ifft(T.*psi_k,[],1),[],2),1),2);
%         psi(:,:,i+1) = ifft(ifft(T_prop.*G_prop2(:,:,i).*fft(fft(G_prop1(:,:,i).*V_prop(:,:,i).*...
%             ifft(ifft(T_prop.*psi_k,[],1),[],2),[],1),[],2),[],1),[],2);
%         T_expect = sum(sum(conj(psi(:,:,i)).*ifft(ifft(T.*psi_k,[],1),[],2),1),2);
        
        V_expect = sum(sum(conj(psi(:,:,i)).*V(:,:,i).*psi(:,:,i),1),2);
        E(i) = real(T_expect + V_expect);
        
        if ~mod(i,10)
            if animate==1
                figure(1);
                clf;
                subplot(2,1,2);
                hold on;
                box on;
                plot(z*1e6,sum(abs(psi(:,:,i)).^2,1)*10)
                plot(z*1e6,sum(V(:,:,i),1)/depth(1,1)/Nx);
                hold off;
                ylim([-2 2])
                xlabel('z (um)');
                subplot(2,1,1);
                hold on;
                box on;
                plot(x*1e6,sum(abs(psi(:,:,i)).^2,2)*10)
                plot(x*1e6,sum(V(:,:,i),2)/depth(1,1)/Nz);
                hold off;
                xlabel('x (um)');
                ylim([-2 2])
                drawnow();
                
            elseif animate==2
                zdata = abs(psi(:,:,i)).^2;
                zdata2 = V(:,:,i)/depth(j,j);
                figure(2);
                clf;
                hold on;
                box on;
                contourf(z*1e6,x*1e6,zdata2,10);
                im = imagesc(z*1e6,x*1e6,zdata);
                im.AlphaData = zdata/max(abs(zdata(:)));
                caxis([-1.2 0])
                hold off;
                shading flat
                xlabel('axial (um)')
                ylabel('radial (um)')
                axis image
                colormap gray
                drawnow();
            end
            
            if save_gif
                h = gcf;
                axis tight manual % this ensures that getframe() returns a consistent size
                %             for n = 1:0.5:5
                % Draw plot for y = x.^n
                %                 x = 0:0.01:1;
                %                 y = x.^n;
                %                 plot(x,y)
                %                 drawnow
                % Capture the plot as an image
                frame = getframe(h);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                if first_save
                    filename = ['merge' num2str(j) '_' datestr(now,'YYmmDD_HHMMSS') '.gif'];
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.05);
                    first_save = false;
                    disp(filename);
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.05);
                end
            end
            
        end
    end
    
%     xcom = reshape(sum(x.*sum(abs(psi).^2,2),1)./sum(sum(abs(psi).^2,2),1),size(t));
%     zcom = reshape(sum(z'.*sum(abs(psi).^2,1),2)./sum(sum(abs(psi).^2,1),2),size(t));
    
    E = (E(:)+depth(ind{j,2}(1),ind{j,2}(2)))/hbar/(omega_x(j,2)/2 + omega_z(j,2)/2);
    E_out(j) = E(end);
    gs_overlap(j) = abs(sum(sum(conj(psi_gs(:,:,2)).*psi(:,:,end),1),2)).^2;
end

end

function out = minjerk(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<=0) = x0;
out(trel>=1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0)*(10*trel(ind).^3 - 15*trel(ind).^4 + 6*trel(ind).^5);
end

function out = minjerk_deriv(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
ind = trel>=0 & trel<=1;
out(ind) = (x1-x0)/(t1-t0)*(30*trel(ind).^2 - 60*trel(ind).^3 + 30*trel(ind).^4);
end

function out = minjerk_accel(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
ind = trel>=0 & trel<=1;
out(ind) = (x1-x0)/(t1-t0)^2*(60*trel(ind).^1 - 180*trel(ind).^2 + 120*trel(ind).^3);
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

