function x = adaptive_grid(W,m,xrange,Nx,radial_boo)

if nargin<5
    radial_boo = 1;
end

xmin = xrange(1);
xmax = xrange(2);

%% set up envelope potential
Wtest = W(xmin);
Nchn = size(Wtest,1);
Wmindiag1 =@(x) reshape(min(diag_nd(W(reshape(x,1,1,numel(x)))),[],1),size(x));
Wmaxdiag1 =@(x) reshape(max(diag_nd(W(reshape(x,1,1,numel(x)))),[],1),size(x));
for i = 1:Nchn
    xe_chn(i) = fminbnd(@(x) mat_el_nd(W(x),i,i),xmin,xmax);
    Emin_chn(i) = mat_el_nd(W(xe_chn(i)),i,i);
end
Emin = Emin_chn(Emin_chn==min(Emin_chn));
xe = xe_chn(Emin_chn==min(Emin_chn));

if radial_boo
    Wenv =@(x) Wmindiag1(x).*(x>=xe) + Emin.*(x<xe);
    Emax = Wmaxdiag1(fminbnd(@(x) -Wmaxdiag1(x),xe,xmax));
else
    Wenv =@(x) Wmindiag1(x);
    Emax = Wmaxdiag1(fminbnd(@(x) -Wmaxdiag1(x),xmin,xmax));
end
Emax = Emax + 1e-6*(Emax-Emin);

%% define adaptive grid
if radial_boo
    xtemp = logspace(log10(xmin),log10(xmax),1e4);
else
    xtemp = linspace(xmin,xmax,1e4);
end

psq_loc_max =@(x) 2*m*(Emax-Wenv(x));
pmax = sqrt(2*m*(Emax-Emin));

% grid mapping jacobian
Jfun =@(x) real(pmax./sqrt(psq_loc_max(x)));

% q is the mapped grid, just a dummy variable here really -- could probably
% eliminate entirely
qmax = integral(@(R) 1./Jfun(R),xmin,xmax);
q = linspace(0,qmax,Nx);

% use interpolation to map the uniformly spaced variable x back to the
% (now unequally spaced) physical variable R
Jtemp = Jfun(xtemp);
qtemp = cumtrapz(xtemp,1./Jtemp);
x = interp1(qtemp,xtemp,q,'pchip','extrap');

end