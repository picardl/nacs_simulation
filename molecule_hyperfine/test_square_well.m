function out = test_square_well()

const = constants();

%% simulation parameters
Nx = 1e3;
rmin = 0;
rmax = 1;
Erange = [0 10];

W =@(x) (x<0.1)*100 + (x>0.9)*100;

%% call the solver
[E_out,nodes_out,psi,r] = cc_logderiv_adaptive_multi([rmin rmax],Nx,W,Erange,1,0);


end