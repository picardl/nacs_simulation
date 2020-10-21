function [psi,E] = fourier_grid_uniform(V,x,m)

if nargin<3
    m = 1;
end

% kinetic energy operator
D2 = fourier_deriv_matrix(x,2);
T = -D2/(2*m);

% build total Hamiltonian
Nx = numel(x);
Vx = V(x);
H = T + spdiags(Vx(:),1,Nx,Nx);
H = (H+H')/2; % enforce hermiticity

% solve eigenvalue problem
[psi,E] = eig(H);
E = diag(E);

end