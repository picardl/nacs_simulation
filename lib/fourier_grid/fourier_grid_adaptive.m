function [R,psi,E] = fourier_grid_adaptive(V,r_range,Emax,m,Nx)

if nargin<4
    Nx = 1e3;
end
if nargin<3
    m = 1;
end

rmin = r_range(1);
rmax = r_range(2);
re = fminbnd(V,rmin,rmax);

Emin = V(re);
pmax = sqrt(2*m*(Emax-Emin));

% define equally spaced grid via local deBroglie wavelength
Jfun =@(r) pmax./sqrt(2*m*(Emax-V(r)));
xmax = integral(@(R) 1./Jfun(R),rmin,rmax);
x = linspace(0,xmax,Nx);

% use interpolation to map the uniformly spaced variable x back to the
% (now unequally spaced) physical variable R
rtemp = linspace(rmin,rmax,1e6);
Jtemp = Jfun(rtemp);
xtemp = cumtrapz(rtemp,1./Jtemp);
R = interp1(xtemp,rtemp,x,'pchip','extrap');
V_R = V(R);
J = Jfun(R);
clear('rtemp');

% kinetic energy operator
D = fourier_deriv_matrix(x,1);
Jinv12 = spdiags(J(:).^(-0.5),0,Nx,Nx);
Jinv = spdiags(1./J(:),0,Nx,Nx);
T = (-1/(2*m))*Jinv12*D*Jinv*D*Jinv12;

% build total Hamiltonian
H = T + spdiags(V_R(:),1,Nx,Nx);
H = (H+H')/2; % enforce hermiticity

% solve eigenvalue problem
[psi,E] = eig(H);
psi = real(psi);
E = diag(E);

% find and delete ghost states
ghost = abs(psi(end,:)).^2 > 0.01;
nghost = sum(ghost);
E(ghost) = [];
psi(:,ghost) = [];
fprintf('%d ghost(s) found and deleted. \n',nghost)

end