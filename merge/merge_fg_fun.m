function [psi,E,Plost,H] = merge_fg_fun(t_move_ramp,power,separation,waist)

if nargin<4
    waist(1) = 0.623e-6;
    waist(2) = 1.064e-6;
end
if nargin<3
    separation = 5e-6;
end
if nargin<2
    power = [1.8e-3 6.8e-3];
end
if nargin<1
    t_move_ramp = [0.04e-3 0.04e-3];
end

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

x0 = 0;
x1 = separation;

xmin = x0 - 3*waist(1);
xmax = x1 + 3*waist(2);
Nx = 400;
x = linspace(xmin,xmax,Nx)';

t0 = 0.0;
t1 = t_move_ramp(1);
t2 = t1 + t_move_ramp(2);

Nt = 400;
t = linspace(0,1.5*t2,Nt);
dt = t(2)-t(1);

%% define hamiltonian
gaussian =@(x0,waist,x) exp(-(x-x0).^2/waist.^2);
d2_dx2 = fourier_deriv_matrix(x,2); % fourier-based derivative matrix

for i = 1:2
    potential = - depth(i,1) * linramp(1,0,t1,t2,t).*(gaussian( minjerk(x0,x1,t0,t1,t),waist(1),x ))...
        - depth(i,2) * gaussian(x1,waist(2),x);
    potential_mat = reshape(potential,[Nx,1,Nt]).*eye(Nx);
    
    kinetic_mat = -(hbar^2/(2*m(i)))*d2_dx2; % kinetic energy operator
    
    H = kinetic_mat + potential_mat;
    
    %% go to basis of initial eigenvectors
    [vecs_init,~] = eigenshuffle(H(:,:,1));
    vecs_init = fliplr(vecs_init);
    
    [vecs_fin,vals_fin] = eigenshuffle(H(:,:,end));
    vecs_fin = fliplr(vecs_fin);
    vals_fin = flipud(vals_fin);
    
    if contains(version,'R2020b')
        H = pagemtimes(pagemtimes(vecs_init',H),vecs_init);
    else
        H = mtimesx(mtimesx(vecs_init',H),vecs_init); % go to basis of initial eigenstates
    end
    
    H = (H + conj(permute(H,[2 1 3])))/2; % enforce hermiticity
    
    %%
    prop = zeros(Nx,Nx,Nt-1);
    parfor j = 1:Nt-1
        prop(:,:,j) = expm(-1i*H(:,:,j)*dt/hbar);
    end
    
    psi = zeros(Nx,1,Nt);
    psi(1,1) = 1;
    for j = 1:Nt-1
        psi(:,1,j+1) = prop(:,:,j)*psi(:,1,j);
    end
    
    psi = psi./sum(abs(psi).^2,1);
    
    psi_finbasis{i} = mtimesx(vecs_fin'*vecs_init,psi);
    
    E{i} = mtimesx(psi,'c',mtimesx(H,psi));
    
    Plost(i) = sum(abs(psi_finbasis(vals_fin>0,1,end)).^2);
    
    figure(1);
    for j = 1:Nt
        clf;
        hold on;
        plot(x*1e6,abs(vecs_init(:,:,1)*psi(:,1,j)).^2)
        plot(x*1e6,potential(:,j)/abs(depth(i,i)))
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

function out = linramp(x0,x1,t0,t1,t)
trel = (t - t0)/(t1-t0);
out = 0*t;
out(trel<0) = x0;
out(trel>1) = x1;
ind = trel>=0 & trel<=1;
out(ind) = x0 + (x1-x0).*trel(ind);
end