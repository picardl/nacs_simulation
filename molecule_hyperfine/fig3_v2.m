clear;

pc = plot_colors();

% figure(3);
% clf;

% data = load('data/fig3b.mat');
% x = data.x*1e6;
% y = data.y;
% dy = data.dy;
% ymin = y(1);
% dymin = dy(1);
% x(1) = [];
% y(1) = [];
% dy(1) = [];
% 
% data = load('data/fig3a.mat');
% x2off = 110.1;
% x2 = data.x - x2off;
% y2 = data.y;
% dy2 = data.dy;
% t = 2.9;

data = load('data/fig3b_alt.mat');
x = data.x(:,3)*1e6;
y = data.y(:,3);
dy = data.dy(:,3);
ymin = 0.35;
dymin = 0.03;

data = load('data/fig3a_alt.mat');
x2off = 249.8; % MHz, freq offset
x2 = data.x(1:end-1,1) - x2off;
y2 = data.y(1:end-1,3);
dy2 = data.dy(1:end-1,3);
t = 2.1;

rabiLineDephase =@(Gamma,gamma,Omega,delta,t) exp(-Gamma*t).*(1 - 0.5*Omega.^2./(Omega.^2 + delta.^2).*(1 - exp(-gamma*t).*cos(sqrt(Omega.^2+delta.^2)*t)));

b = [   
    y(1)-ymin
    ymin
    0.03
    0.03
    0.19
    0
    ]';

rabi_vs_t =@(b,t) b(1) + b(2)*rabiLineDephase(b(3),b(4),2*pi*b(5),0,t);
rabi_vs_freq =@(b,delta) b(1) + b(2)*rabiLineDephase(b(3),b(4),2*pi*b(5),2*pi*(delta-b(6)),t);

F =@(b) cat(1,(b(1)-ymin)./dymin,(rabi_vs_t(b,x) - y)./dy,(rabi_vs_freq(b,x2) - y2)./dy2);

[b,~,resid,~,~,~,J] = lsqnonlin(F,b);
CI = nlparci(b,resid,'jacobian',J,'alpha',0.6827);
db = diff(CI,[],2)'/2;

Patom = b(1);
Pmol = b(2);
Gamma = b(3);
gamma = b(4);
rabi = b(5);
delta0 = b(6);

Patom_err = db(1);
Pmol_err = db(2);
Gamma_err = db(3);
gamma_err = db(4);
rabi_err = db(5);
delta0_err = db(6);

xfit = linspace(min(x),28,5e2)';
yfit = rabi_vs_t(b,xfit);

% pi time
[tpi,tpi_err] = prop_errors(@(x) 1/(2*x),rabi,rabi_err);

% minimum point of fb mol pop on first rabi flop
[ym,dym] = nlpredci(rabi_vs_t,tpi,b,resid,'jacobian',J);

% maximum point at 2pi time
[yreturn,dyreturn] = nlpredci(rabi_vs_t,2*tpi,b,resid,'jacobian',J);

% 1-way efficiency is total remaining population patom+pmol*exp(-gamma*tpi)
% minus minimum point of rabi flop
[eff_1way,eff_1way_err] = prop_errors(@(patom,pmol,Gamma,tpi,ybot) (patom + pmol*exp(-Gamma*tpi))-ybot,[Patom Pmol Gamma tpi ym],[Patom_err Pmol_err Gamma_err tpi_err dym]);

% population lost by the 2pi time
[lost,lost_err] = prop_errors(@(Gamma,tpi,pmol) pmol*(1-exp(-Gamma*2*tpi)),[Pmol Gamma tpi],[Pmol_err Gamma_err tpi_err]);

% population remaining in g.s. due to decoherence at 2pi time. Slightly
% non-obvious expression, derived in mathematica "raman_rabi_efficiency.nb"
[remain,remain_err] = prop_errors(@(pmol,Gamma,gamma,Omega) pmol*0.5*exp(-2*pi*(gamma+Gamma)./Omega).*(exp(2*pi*gamma./Omega)-1),[Pmol Gamma gamma 2*pi*rabi],[Pmol_err Gamma_err gamma_err 2*pi*rabi_err]);

% round-trip efficiency: population that went to the g.s. and back
[eff_rt,eff_rt_err] = prop_errors(@(pmol,remain,lost) pmol-remain-lost,[Pmol,remain,lost],[Pmol_err,remain_err,lost_err]);

[Omega_Gamma_ratio,Omega_Gamma_ratio_err] = prop_errors(@(Omega,Gamma) 2*pi*Omega/Gamma,[rabi Gamma],[rabi_err Gamma_err]);

s = errstr(b,db,2);
fprintf('ymin = %s\ncontrast = %s\nGamma = %s\ngamma = %s\nRabi = %s\nDelta0 = %s\n',s{:})

[transfer_frac,transfer_frac_err] = prop_errors(@(a,b) a./b,[eff_1way b(2)],[eff_1way_err db(2)]);
[fb_frac,fb_frac_err] = prop_errors(@(a,b,c) (a-b),[y(1) ymin b(2)],[dy(1) dymin db(2)]);

fprintf('transfer fraction = %s\n',errstr(transfer_frac,transfer_frac_err))

fprintf('1 way total = %s\n',errstr(eff_1way,eff_1way_err))
fprintf('round trip total = %s\n',errstr(eff_rt,eff_rt_err))
fprintf('Omega/Gamma = %s\n',errstr(Omega_Gamma_ratio,Omega_Gamma_ratio_err))

% [g,~,~,cov] = nlinfit(x2,y2,rabi_vs_freq,b,'Weights',dy2.^-2);
% dg = reshape(sqrt(diag(cov)),[3 1]);
xfit2 = linspace(-1.5,1.5,3e2);
yfit2 = rabi_vs_freq(b,xfit2);

figure(3);
clf;
set(gcf,'units','inches','position',[5 5 3.4 2.2]);
hold on; box on;

xp = [-2 30];
curve1 = [1 1]*(ymin+dymin);
curve2 = [1 1]*(ymin-dymin);
plot(xp, curve1,'color','none');
plot(xp, curve2,'color','none');
xp2 = [xp, fliplr(xp)];
inBetween = [curve1, fliplr(curve2)];
fill1 = fill(xp2, inBetween, 'g');
set(fill1','facecolor',[1 1 1]*0.8,'edgecolor','none')
plot(xp,[1 1]*ymin,'--k');

xp2 = [-2 30];
curve3 = [1 1]*(y(1)+dy(1));
curve4 = [1 1]*(y(1)-dy(1));
plot(xp2, curve3,'color','none');
plot(xp2, curve4,'color','none');
x3 = [xp2, fliplr(xp2)];
inBetween = [curve3, fliplr(curve4)];
fill2 = fill(x3, inBetween, 'g');
set(fill2','facecolor',[1 1 1]*0.8,'edgecolor','none')
plot(xp2,[1 1]*y(1),'--k');

plot(xfit,yfit,'-k','linewidth',0.75);
set(errorbar(x,y,dy,'o'),'color',pc.lc{2},'markersize',4,'markerfacecolor',pc.fc{2},'markeredgecolor',pc.ec{2},'capsize',2,'linewidth',0.75);
set(errorbar(0,ymin,dymin,'v'),'color',pc.lc{1},'markersize',4,'markerfacecolor',pc.fc{1},'markeredgecolor',pc.ec{1},'capsize',2,'linewidth',0.75)
hold off;

% xlim([-1 26])
xlim([-0.5 8.5])
ylim([0.29 0.85])
% ylim([0.28 0.75])

xlabel('Pulse duration (µs)')
ylabel('Na + Cs survival')

set(gca,'fontsize',8,'color','none','linewidth',0.75,'xminortick','on','yminortick','on');

ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom+0.01*ax_height ax_width*0.98 ax_height*1.0];
ax.Position = [left bottom+0.01*ax_height ax_width*0.98 ax_height*0.98];

% set(annotation('textbox', [0.10, 0.4, 0, 0], 'string', '(b)'),'fontsize',10)

% save
fname = 'fig3c.eps';
saveas(gcf,['../figures/' fname],'epsc');
fprintf('saved file %s\n',fname);

%%

% data = load('data/fig3a.mat');
% x2 = data.x;
% y2 = data.y;
% dy2 = data.dy;

% rabi_fun=@(b,x) b(1) + b(2)*(1-rabiLine(2*pi*(x-b(3))*1e6,t,2*pi*b(4)*1e6));

% guess = [0.3 0.35 0 0.3];
% [b,~,~,cov] = nlinfit(x2,y2,rabi_fun,guess,'Weights',dy2.^-2);
% db = reshape(sqrt(diag(cov)),size(b));

figure(4);
clf;
% set(gcf,'units','inches','position',[5 5 3.4/2 1.8]);
set(gcf,'units','inches','position',[5 7.5 3.4 2.2]);

% yfit2 = rabi_solve([b(1:2) 2*pi*0 0],2*pi*(xfit-b(3)*100)*1e6,2*pi*b(4)*1e6,2.7e-6);

% ax = axes('units','normalized','position',[0.55 0.62 0.4 0.32]);
hold on;
box on;
xp = [-2 30];
curve1 = [1 1]*(ymin+dymin);
curve2 = [1 1]*(ymin-dymin);
plot(xp, curve1,'color','none');
plot(xp, curve2,'color','none');
xp2 = [xp, fliplr(xp)];
inBetween = [curve1, fliplr(curve2)];
fill1 = fill(xp2, inBetween, 'g');
set(fill1','facecolor',[1 1 1]*0.8,'edgecolor','none')
plot(xp,[1 1]*ymin,'--k');

xp2 = [-2 10];
curve3 = [1 1]*(y(1)+dy(1));
curve4 = [1 1]*(y(1)-dy(1));
plot(xp2, curve3,'color','none');
plot(xp2, curve4,'color','none');
x3 = [xp2, fliplr(xp2)];
inBetween = [curve3, fliplr(curve4)];
fill2 = fill(x3, inBetween, 'g');
set(fill2','facecolor',[1 1 1]*0.8,'edgecolor','none')
plot(xp2,[1 1]*y(1),'--k');

plot(xfit2,yfit2,'-k','linewidth',0.75);
set(errorbar(x2,y2,dy2,'o'),'color',pc.lc{2},'markersize',4,'markerfacecolor',pc.fc{2},'markeredgecolor',pc.ec{2},'capsize',2,'linewidth',0.75)
hold off;
%
xlim([-1 1])
ylim([0.29 0.85])
% xlim([-1 1])

% ylim([0.28 0.75])
%
xlabel('Stokes detuning (MHz)','fontsize',8)
ylabel('Na + Cs survival')
%
ax = gca;
set(ax,'fontsize',8,'color','none','linewidth',0.75,'xminortick','on','yminortick','on');

ax.Position = [left bottom+0.01*ax_height ax_width*0.98 ax_height*0.98];
%
% set(gcf,'units','inches','position',[5 5 2.6 2]);
% ax.Position = [left+0.01*ax_width bottom+0.01*ax_height ax_width*0.98 ax_height*0.98];

% set(annotation('textbox', [0.10, 0.4, 0, 0], 'string', '(a)'),'fontsize',10)

%% save
% fname = 'fig3b.eps';
% saveas(gcf,['../figures/' fname],'epsc');
% fprintf('saved file %s\n',fname);


