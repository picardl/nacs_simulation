
function [E_out,nodes_out,psi,x] = cc_logderiv_adaptive_multi(xrange,Nx,W,Erange,m,radial_boo,prop_wfn,verbose,plot_boo)

if nargin<9
    plot_boo = 0;
end
if nargin<8
    verbose = 0;
end
if nargin<7
    prop_wfn = 1;
end
if nargin<6
    radial_boo = 1;
end

Wtest = W(xrange(1));
Nchn = size(Wtest,1);

x = adaptive_grid(W,m,xrange,Nx,radial_boo);

%% upwards and downwards propagation initialization
xc = (x(2:end)+x(1:end-1))/2;
xabc = zeros(1,2*Nx-1);
xabc(1:2:end) = x;
xabc(2:2:end) = xc;
h = diff(x)/2;
Wabc = W(reshape(xabc,1,1,[]));

%% diagonalize in each sector
[Vabc,D] = eigenshuffle(Wabc(:,:,1:2:end));
Wb = eye(Nchn).*reshape(D,[1,Nchn,Nx]);
Wac = Wabc(:,:,2:2:end-1);
Wac = pagemtimes(conj(permute(Vabc(:,:,1:end-1),[2,1,3])),pagemtimes(Wac,Vabc(:,:,1:end-1)));
Wb_copy = Wabc(:,:,1:2:end);
Wabc = 0*Wabc;
Wabc(:,:,2:2:end-1) = Wac;
Wabc(:,:,1:2:end) = Wb;

if plot_boo
    Wplot = plot_W();
end

%% plot node counting and Ymatch
% 
% Nsamp = 400;
% Esamp = linspace(Erange(1),Erange(2),Nsamp);
% 
% nodes = zeros(1,Nsamp);
% evals = zeros(Nchn,Nsamp);
% 
% for i = 1:Nsamp
%     [~,~,~,nodes(i),evals(:,i)] = get_Ymatch_eval(Esamp(i));
% end
% 
% const = constants();
% 
% figure(101);
% clf;
% subplot(2,1,1);
% plot(Esamp*const.hartree/const.h*1e-12,evals,'.-');
% ylabel('eig(Ymatch)')
% 
% subplot(2,1,2);
% plot(Esamp*const.hartree/const.h*1e-12,nodes,'.-');
% xlabel('E (THz)')
% ylabel('nodes')
% 
% Esamp;

%%
    function [Ymatch,Mu,Md,nodes,evals] = get_Ymatch_eval(E)
        
        kup = sqrt(2*m*(diag(Wabc(:,:,1))-E));
        Yiup = diag(abs(kup));
        kdn = sqrt(2*m*(diag(Wabc(:,:,end))-E));
        Yidn = -diag(abs(kdn));
        
        node_boo = zeros(1,Nx);
        
        Y = zeros(Nchn,Nchn,Nx+1);
        Y(:,:,1) = Yiup;
        Y(:,:,end) = Yidn;
        M = zeros(Nchn,Nchn,Nx+1);
        
        Ediff = E-diag_nd(Wabc(:,:,1:2:end));
        chns_allowed = any(Ediff > 0,2);
        
        j = 1;
        stop = false;
        while ~stop
            [Y(:,:,Nx+1-j),M(:,:,Nx+1-j),~] = modlogderiv_propagator_v4(Y(:,:,Nx-j+2),...
                Wabc(:,:,2*Nx - 2*j + (-1:1)),Vabc(:,:,Nx-j),h(Nx-j),E,m,1);
            if (any(eig(Y(:,:,Nx+1-j))>0) && all((Ediff(:,Nx+1-j)>0) == chns_allowed)) || j==Nx-1
                stop = true;
            end
            j = j+1;
        end
        Nd = j;
        Nu = Nx-Nd+1;
        
        for j = 1:(Nu-1)
            [Y(:,:,j+1),M(:,:,j+1),node_boo(j)] = modlogderiv_propagator_v4(Y(:,:,j),...
                Wabc(:,:,2*j + (-1:1)),Vabc(:,:,j),h(j),E,m,0);
        end
        Ycopy = Y;
        
        for j = Nu:(Nx-1)
            [Ycopy(:,:,j+1),~,node_boo(j)] = modlogderiv_propagator_v4(Ycopy(:,:,j),...
                Wabc(:,:,2*j + (-1:1)),Vabc(:,:,j),h(j),E,m,0);
        end
        
        Ymatch = Y(:,:,Nu+1)-Y(:,:,Nu);
        Mu = M(:,:,1:Nu);
        Md = M(:,:,(Nu+1):end);
        nodes = sum(node_boo);
        
        [evecs,evals] = eig(Ymatch);
        evals = real(diag(evals));
        [~,ind] = max(abs(evecs).^2,[],1);
        [~,ord] = sort(ind);
        evals = evals(ord);
    end

%% use bisection to find regions with only one eigenvalue
if verbose
    fprintf('Finding energy ranges with single eigenvalues...')
end

[Ymatch(:,:,1),~,~,nodes_samp(1),Ymatch_evals(:,1)] = get_Ymatch_eval(Erange(1));
[Ymatch(:,:,2),~,~,nodes_samp(2),Ymatch_evals(:,2)] = get_Ymatch_eval(Erange(2));

unodes = nodes_samp(1):nodes_samp(2);

i = 1;
Elist = Erange;
while ~all(ismember(unodes,nodes_samp))
    if (nodes_samp(i+1)-nodes_samp(i))>1
        Emid = (Elist(i+1)+Elist(i))/2;
        [Ym,~,~,nodes_mid,evals_mid] = get_Ymatch_eval(Emid);
        
        Elist = [Elist(1:i) Emid Elist(i+1:end)];
        nodes_samp = [nodes_samp(1:i) nodes_mid nodes_samp(i+1:end)];
        Ymatch = cat(3,Ymatch(:,:,1:i),Ym,Ymatch(:,:,i+1:end));
        Ymatch_evals = cat(2,Ymatch_evals(:,1:i),evals_mid,Ymatch_evals(:,i+1:end));
    else
        i = i+1;
    end
    
    if plot_boo
        plot_Elist()
    end
end

if verbose
    fprintf('done.\n')
end

%% find the sign changes of the objective function
if verbose
    fprintf('Finding sign changes of objective function...\n')
end

itermax = 100;
tol = 1e-4;
E_out = [];
nodes_out = [];
evals_out = [];
psi = [];
for i = 1:numel(unodes)-1
    
    if verbose
        fprintf('\n')
        fprintf('step %d of %d\n',i,numel(unodes)-1)
        fprintf('nodes: %d\n',unodes(i))
    end
    
    iter = 1;
    while true
        kk = find(nodes_samp==unodes(i),1,'last');
        Elower = Elist(kk);
        Eupper = Elist(kk+1);
        Emid = (Eupper + Elower)/2;
        
        [Ym,Mu,Md,nodes_mid,evals_mid] = get_Ymatch_eval(Emid);
        
        if (min(abs(evals_mid)) < tol)
            E_out = cat(2,E_out,Emid);
            nodes_out = cat(2,nodes_out,nodes_mid);
            evals_out = cat(2,evals_out,evals_mid);
            
            Nd = size(Md,3);
            Nu = size(Mu,3);
            
            %% propagate the wavefunction
            if prop_wfn
                if verbose
                    fprintf('propagating wfn...\n')
                end
                
                [V,D] = eig(Ym);
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
                
                psi_new = real(cat(2,psi_up,psi_dn(:,2:end)));
                psi_new = psi_new/sqrt(trapz(x,sum(abs(psi_new).^2,1)));
                
                psi = cat(3,psi,psi_new);
                
                if plot_boo
                    plot_psi();
                end
            end
            
            break;
        elseif iter>itermax
            print('failed to converge on eval %d\n',unodes(i))
            break;
        end
        
        
        Elist = [Elist(1:kk) Emid Elist(kk+1:end)];
        nodes_samp = [nodes_samp(1:kk) nodes_mid nodes_samp(kk+1:end)];
        Ymatch_evals = cat(2,Ymatch_evals(:,1:kk),evals_mid,Ymatch_evals(:,kk+1:end));
        
        if plot_boo
            plot_Elist()
        end
        iter = iter + 1;
    end
end

if verbose
    fprintf('done.\n')
end

%% plot functions
    function plot_Elist()
        figure(2); clf;
        
        sp(1) = subplot(2,1,1);
        plot(Elist,Ymatch_evals','.-');
        ylim([-100 100])
        xlabel('E (E_h)');
        ylabel('eig(Y_{match})')
        xlim(Erange)
        
        sp(2) = subplot(2,1,2);
        plot(Elist,nodes_samp,'.-k');
        ylabel('nodes')
        xlabel('E (E_h)');
        
        xlim(Erange)
        linkaxes(sp,'x');
        
        
        drawnow();
    end

    function plot_psi()
        psi_plot = 0.5*(psi./max(max(abs(psi),[],1),[],2))*peak2peak(Erange)/numel(unodes) + reshape(E_out,1,1,[]);
        
        figure(1);
        clf;
        Wplot = plot_W();
        hold on;
        for psi_plot_ind = 1:size(psi,3)
            p = plot(x,psi_plot(:,:,psi_plot_ind));
            arrayfun(@(p,w) set(p,'Color',get(w,'Color')),p,Wplot)
        end
        hold off;
        drawnow();
    end

    function Wplot = plot_W()
        if Nchn==1
            W_adiab = Wb(:)';
            W_diab = Wb_copy(:)';
        else
            W_adiab = diag_nd(Wb);
            W_diab = diag_nd(Wb_copy);
        end
        
        figure(1);
        clf;
        hold on;
        box on;
        Wplot = plot(x,W_adiab);
        diabplot = plot(x,W_diab,'--');
        arrayfun(@(a,b) set(a,'color',get(b,'color')),diabplot,Wplot);
        plot(xrange,[1 1]*Erange(1),'-k')
        plot(xrange,[1 1]*Erange(2),'-k')
        hold off;
        set(gca,'xscale','log')
        xlabel('R (a_0)')
        ylabel('E (Hartree)')
        xlim(xrange)
        ylim(Erange + [-1 1]*peak2peak(Erange)*0.25)
    end

end

