
function [E_out,nodes_out,err_est,psi,x] = cc_logderiv_adaptive_multi(xrange,Nx,W,Erange,m,radial_boo,prop_wfn,verbose,plot_boo,richardson_extrap)

if nargin<10
    richardson_extrap = 1; % currently semi-broken
end
if nargin<9
    plot_boo = 1;
end
if nargin<8
    verbose = 1;
end
if nargin<7
    prop_wfn = 1;
end
if nargin<6
    radial_boo = 1;
end

Escale = abs(mean(Erange));
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
% Emin = Emin(1);
xe = xe_chn(Emin_chn==min(Emin_chn));
% xe = xe(1);

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

%% upwards and downwards propagation initialization
xc = (x(2:end)+x(1:end-1))/2;
xabc = zeros(1,2*Nx-1);
xabc(1:2:end) = x;
xabc(2:2:end) = xc;
h = diff(x)/2;
Wabc = W(reshape(xabc,1,1,[]));

%% dense sampling for richardson extrapolation
Nx_dense = numel(xabc);
x_dense = xabc;
xc_dense = (x_dense(2:end)+x_dense(1:end-1))/2;
xabc_dense = zeros(1,2*Nx_dense-1);
xabc_dense(1:2:end) = x_dense;
xabc_dense(2:2:end) = xc_dense;
h_dense = diff(x_dense)/2;
Wabc_dense = zeros(Nchn,Nchn,2*Nx_dense-1);
Wabc_dense(:,:,1:2:end) = Wabc;
Wabc_dense(:,:,2:2:end) = W(reshape(xc_dense,1,1,[]));

Nx_dense2 = numel(xabc_dense);
x_dense2 = xabc_dense;
h_dense2 = diff(x_dense2)/2;
xc_dense2 = (x_dense2(2:end)+x_dense2(1:end-1))/2;
Wabc_dense2 = zeros(Nchn,Nchn,2*Nx_dense2-1);
Wabc_dense2(:,:,1:2:end) = Wabc_dense;
Wabc_dense2(:,:,2:2:end) = W(reshape(xc_dense2,1,1,[]));

%%
fun_evals_tot = 0;

Ne = 1e2;
Esamp = linspace(Erange(1),Erange(2),Ne);

Ymatch = zeros(Nchn,Nchn,Ne);
nodes = zeros(Nchn,Ne);
evals = zeros(Nchn,Ne);
for i = 1:Ne
    [~,Ymatch(:,:,i),~,~,nodes(:,i),evals(:,i)] = get_Ymatch_eval(Esamp(i),Wabc,Nx,h);
end

figure(1001);
clf;
subplot(2,1,1);
plot(Esamp,nodes');

subplot(2,1,2);
plot(Esamp,evals);
ylim([-100 100])

Esamp;


%%
    function [out,Ymatch,Mu,Md,nodes,evals] = get_Ymatch_eval(E,Wabc,Nx,h,eval_ind)
        
        kup = sqrt(2*m*(diag(Wabc(:,:,1))-E));
        Yiup = diag(abs(kup));
        kdn = sqrt(2*m*(diag(Wabc(:,:,end))-E));
        Yidn = -diag(abs(kdn));
        
        node_boo = zeros(Nchn,Nx);
        
        Y = zeros(Nchn,Nchn,Nx+1);
        Y(:,:,1) = Yiup;
        Y(:,:,end) = Yidn;
        M = zeros(Nchn,Nchn,Nx+1);
        
        Ediff = E-diag_nd(Wabc(:,:,1:2:end));
        chns_allowed = any(Ediff>0,2);
        
        
        j = 1;
        stop = false;
        while ~stop
            [Y(:,:,Nx+1-j),M(:,:,Nx+1-j),node_boo(:,Nx+1-j)] = modlogderiv_propagator_v3(Y(:,:,Nx-j+2),...
                Wabc(:,:,2*Nx - 2*j + (-1:1)),h(Nx-j),E,m,1);
            if (any(eig(Y(:,:,Nx+1-j))>0) && all((Ediff(:,Nx+1-j)>0) == chns_allowed)) || j==Nx-1
                stop = true;
            end
            j = j+1;
        end
        Nd = j;
        Nu = Nx-Nd+1;
        
        for j = 1:(Nu-1)
            [Y(:,:,j+1),M(:,:,j+1),node_boo(:,j)] = modlogderiv_propagator_v3(Y(:,:,j),...
                Wabc(:,:,2*j + (-1:1)),h(j),E,m,0);
        end
        
        Ymatch = Y(:,:,Nu+1)-Y(:,:,Nu);
        Mu = M(:,:,1:Nu);
        Md = M(:,:,(Nu+1):end);
        nodes = sum(sum(node_boo,2),1);
        
%         Nchn_allowed = sum(chns_allowed);
        
        
%         any(ev<0) && all((Ediff(:,j)>0) == chns_allowed)
        
%         j = 1;
%         stop = false;
%         while ~stop
%             [Y(:,:,j+1),M(:,:,j+1),node_boo(:,j)] = modlogderiv_propagator_v3(Y(:,:,j),...
%                 Wabc(:,:,2*j + (-1:1)),h(j),E,m,0);
%             ev = eig(Y(:,:,j+1));
%             j = j+1;
%             if any(ev<0) && all((Ediff(:,j)>0) == chns_allowed)
%                 stop = true;
%             end
%         end
%         Nu = j;
%         Nd = Nx-Nu;
%         
%         for j = 1:Nd
%             [Y(:,:,Nx+1-j),M(:,:,Nx+1-j),node_boo(:,Nx+1-j)] = modlogderiv_propagator_v3(Y(:,:,Nx-j+2),...
%                 Wabc(:,:,2*Nx - 2*j + (-1:1)),h(Nx-j),E,m,1);
%         end
%         node_boo_copy = node_boo;
        
%         for j = (Nd+1):(Nx-1)
%             [~,~,node_boo_copy(:,Nx+1-j)] = modlogderiv_propagator_v3(Y(:,:,Nx-j+2),...
%                 Wabc(:,:,2*Nx - 2*j + (-1:1)),h(Nx-j),E,m,1);
%         end
        
%         Ymatch = Y(:,:,Nu+1)-Y(:,:,Nu);
% %         Ymatch = Y(:,:,Nu+1)+Y(:,:,Nu);
%         Mu = M(:,:,1:Nu);
%         Md = M(:,:,(Nu+1):end);
%         nodes = sum(sum(node_boo,2),1);
        
        [evecs,evals] = eig(Ymatch);
        evals = real(diag(evals));
        [~,ind] = max(abs(evecs).^2,[],1);
        [~,ord] = sort(ind);
        evals = evals(ord);
        if nargin<5
            [~,ind] = min(abs(evals));
            out = evals(ind);
        else
            out = evals(eval_ind);
        end
        fun_evals_tot = fun_evals_tot+1;
    end

%% use bisection to find regions with only one eigenvalue
if verbose
    fprintf('Finding energy ranges with single eigenvalues...')
end

[~,Ymatch(:,:,1),~,~,nodes_samp(1),Ymatch_evals(:,1)] = get_Ymatch_eval(Erange(1),Wabc,Nx,h);
[~,Ymatch(:,:,2),~,~,nodes_samp(2),Ymatch_evals(:,2)] = get_Ymatch_eval(Erange(2),Wabc,Nx,h);

unodes = nodes_samp(1):nodes_samp(2);

i = 1;
while ~all(ismember(unodes,nodes_samp))
    if (nodes_samp(i+1)-nodes_samp(i))>1
        Emid = (Erange(i+1)+Erange(i))/2;
        [~,Ym,~,~,nodes_mid,evals_mid] = get_Ymatch_eval(Emid,Wabc,Nx,h);
        
        Erange = [Erange(1:i) Emid Erange(i+1:end)];
        nodes_samp = [nodes_samp(1:i) nodes_mid nodes_samp(i+1:end)];
        Ymatch = cat(3,Ymatch(:,:,1:i),Ym,Ymatch(:,:,i+1:end));
        Ymatch_evals = cat(2,Ymatch_evals(:,1:i),evals_mid,Ymatch_evals(:,i+1:end));
    else
        i = i+1;
    end
    
    if plot_boo
        figure(101); clf;
        subplot(2,1,1);
        plot(Erange,Ymatch_evals','.-');
        xlabel('E (E_h)');
        ylabel('eig(Y_{match})')
%         set(gca,'fontsize',14);
        subplot(2,1,2);
        plot(Erange,nodes_samp,'.-k');
        ylabel('nodes')
        xlabel('E (E_h)');
%         set(gca,'fontsize',14);
        drawnow();
    end
end

if verbose
    fprintf('done.\n')
end

%% find the divergences of the objective function
if verbose
    fprintf('Finding sign changes of objective function...')
end

itermax = 6;
cut = [];
for i = 1:numel(unodes)
    
    iter = 1;
    while true
        kk = find(nodes_samp==unodes(i),1,'first');
        Elower = Erange(max(kk-1,1));
        Eupper = Erange(kk);
        Emid = (Eupper + Elower)/2;
        
        % midpoint lower node
        [~,Ym,~,~,nodes_mid,evals_mid] = get_Ymatch_eval(Emid,Wabc,Nx,h);
        if nodes_mid==unodes(i)
            Erange = [Erange(1:kk-1) Emid Erange(kk:end)];
            nodes_samp = [nodes_samp(1:kk-1) nodes_mid nodes_samp(kk:end)];
            Ymatch = cat(3,Ymatch(:,:,1:kk-1),Ym,Ymatch(:,:,kk:end));
            Ymatch_evals = cat(2,Ymatch_evals(:,1:kk-1),evals_mid,Ymatch_evals(:,kk:end));
        else
            Erange = [Erange(1:kk-1) Emid Erange(kk:end)];
            nodes_samp = [nodes_samp(1:kk-1) nodes_mid nodes_samp(kk:end)];
            Ymatch = cat(3,Ymatch(:,:,1:kk-1),Ym,Ymatch(:,:,kk:end));
            Ymatch_evals = cat(2,Ymatch_evals(:,1:kk-1),evals_mid,Ymatch_evals(:,kk:end));
        end
        
        
        k = find(nodes_samp==unodes(i),1,'last');
        Eupper2 = Erange(min(k+1,numel(Erange)));
        Elower2 = Erange(k);
        Emid2 = (Eupper2 + Elower2)/2;
        
        % midpoint upper node
        [~,Ym2,~,~,nodes_mid2,evals_mid2] = get_Ymatch_eval(Emid2,Wabc,Nx,h);
        if nodes_mid2==unodes(i)
            Erange = [Erange(1:k) Emid2 Erange(k+1:end)];
            nodes_samp = [nodes_samp(1:k) nodes_mid2 nodes_samp(k+1:end)];
            Ymatch = cat(3,Ymatch(:,:,1:k),Ym2,Ymatch(:,:,k+1:end));
            Ymatch_evals = cat(2,Ymatch_evals(:,1:k),evals_mid2,Ymatch_evals(:,k+1:end));
        else
            Erange = [Erange(1:k) Emid2 Erange(k+1:end)];
            nodes_samp = [nodes_samp(1:k) nodes_mid2 nodes_samp(k+1:end)];
            Ymatch = cat(3,Ymatch(:,:,1:k),Ym2,Ymatch(:,:,k+1:end));
            Ymatch_evals = cat(2,Ymatch_evals(:,1:k),evals_mid2,Ymatch_evals(:,k+1:end));
        end
        
        if plot_boo
            figure(102); clf;
            subplot(2,1,1);
            plot(Erange,Ymatch_evals','.-');
            xlabel('E (E_h)');
            ylabel('eig(Y_{match})')
%             set(gca,'fontsize',14);
            subplot(2,1,2);
            plot(Erange,nodes_samp,'.-k');
            ylabel('nodes')
            xlabel('E (E_h)');
%             set(gca,'fontsize',14);
            drawnow();
        end
        iter = iter + 1;
        if iter > itermax
            if verbose
                fprintf('\nFailed to find a zero crossing for eigenvalue %d. Cutting it :(\n',unodes(i));
            end
            cut = [cut find(nodes_samp==unodes(i))];;

            break;
        end
        
        if size(Ymatch_evals(:,nodes_samp==unodes(i)),2)>2
            if sum(sum(diff(sign(Ymatch_evals(:,nodes_samp==unodes(i))),[],2),2),1)>0
                break;
            end
        end
    end
end

Erange(cut) = [];
nodes_samp(cut) = [];
Ymatch(:,:,cut) = [];
Ymatch_evals(:,cut) = [];
unodes = unique(nodes_samp);

if verbose
    fprintf('done.\n')
end

%% finally, find the zeros using fzero
if verbose
    fprintf('Converging on zeros...')
end

unodes = unique(nodes_samp);
psi_count = 1;
for i = 1:numel(unodes)
    n = nodes_samp==unodes(i);
    E_n = Erange(n);
    
    for k = 1:Nchn
        sign_change_ind = find(diff(sign(Ymatch_evals(k,n)),[],2)~=0);
        if any(sign_change_ind)
            Elo = E_n(sign_change_ind);
            Ehi = E_n(sign_change_ind+1);
            Ezero = Escale*fzero(@(E) get_Ymatch_eval(E*Escale,Wabc,Nx,h,k),[Elo(1) Ehi(1)]/Escale);
            if richardson_extrap
                try
                    Ezero_dense = Escale*fzero(@(E) get_Ymatch_eval(E*Escale,Wabc_dense,Nx_dense,h_dense,k),[Elo Ehi]/Escale);
                    Ezero_dense2 = Escale*fzero(@(E) get_Ymatch_eval(E*Escale,Wabc_dense2,Nx_dense2,h_dense2,k),[Elo Ehi]/Escale);
                    E_extrap1 = (2*Ezero_dense-Ezero);
                    E_extrap2 = (Ezero - 6*Ezero_dense + 8*Ezero_dense2)/3;
                    E_out(psi_count) = E_extrap2;
                    err_est(psi_count) = abs((E_extrap2-E_extrap1)/Nx);
                catch
                    warning('Richarson extrapolation failed at node %d',unodes(i));
                    E_out(psi_count) = Ezero;
                    err_est(psi_count) = 0;
                end
            else
                E_out(psi_count) = Ezero;
                err_est(psi_count) = 0;
            end
            [~,Ymatch_k,Mu,Md,nodes_out(psi_count)] = get_Ymatch_eval(Ezero,Wabc,Nx,h);
            
            Nd = size(Md,3);
            Nu = size(Mu,3);
            
            %% propagate the wavefunction
            if prop_wfn
                [V,D] = eig(Ymatch_k);
                D = diag(D);
                zero_ind = abs(D)==min(abs(D));
                V = V(:,zero_ind);
                psi_up = zeros(Nchn,Nu);
                psi_up(:,end) = V;
                for kk = Nu:-1:2
                    psi_up(:,kk-1) = Mu(:,:,kk)*psi_up(:,kk);
                end
                
                psi_dn = zeros(Nchn,Nd);
                psi_dn(:,1) = V;
                for kk = 1:Nd-1
                    psi_dn(:,kk+1) = Md(:,:,kk)*psi_dn(:,kk);
                end
                
                psi(:,:,psi_count) = real(cat(2,psi_up,psi_dn(:,2:end)));
                psi(:,:,psi_count) = psi(:,:,psi_count)/sqrt(trapz(x,sum(abs(psi(:,:,psi_count)).^2,1)));
                
                if plot_boo
                    figure(103);
                    plot(x,psi(:,:,psi_count),'linewidth',2);
                    set(gca,'xscale','log');
%                     set(gca,'fontsize',14)
                    xlabel('R (a_0)')
                    ylabel('\psi(R)')
                    title(sprintf('E = %g, nodes = %g',E_out(psi_count),nodes_out(psi_count)));
                    drawnow();
                end
                
                psi_count = psi_count+1;
            else
                psi(:,:,psi_count) = [];
            end
        end
    end
    
end

if verbose
    fprintf('done.\n')
end

end

