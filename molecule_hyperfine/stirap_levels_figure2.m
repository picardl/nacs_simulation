function stirap_levels_figure2

c = 299792458;
wavenum2hartree = 4.55634011e-6;
bohr2angstrom = 0.529177210903;

wavenum920 = 325130/c*1e7;
wavenum635 = 472166/c*1e7;

rmin = 4.7;
rmax = 20;
N = 1e3;

Eb = 4954.24;

r = linspace(rmin,rmax,N);

X = NaCsXPES(r)/wavenum2hartree;
a = NaCsaPES(r)/wavenum2hartree;
c = NaCscPES(r)/wavenum2hartree;
B = NaCsBPES(r)/wavenum2hartree;
b = NaCsbbPES(r)/wavenum2hartree;

yoff = min(X(:));
yscale = 1e3;

X = (X - yoff)/yscale;
a = (a - yoff)/yscale;
c = (c - yoff)/yscale;
B = (B - yoff)/yscale;
b = (b - yoff)/yscale;

r = r*bohr2angstrom;

lw = 2;
        
stokes_color = [0.3 0.3 0.7];
pump_color = [0.3 0.3 0.7];

figure(1); clf;
hold on; box off;

% curves
Xplot = plot(r,X,'linewidth',lw);
plot([4.75 9],[1 1]*Eb*1e-3,'-k','linewidth',1)
plot([0 0],[0 0],'color','none')

aplot = plot(r,a,'linewidth',lw);
Bplot = plot(r,B,'linewidth',lw);
bplot = plot(r,b,'linewidth',lw);
cplot = plot(r,c,'linewidth',lw);

% lasers
pump = plot([1 1]*6.5,[Eb wavenum920+Eb]/yscale,'-m','linewidth',2);
stokes = plot([1 1]*3.84,[Eb-wavenum635+wavenum920 wavenum920+Eb]/yscale,'-','linewidth',2);
set(stokes,'color',stokes_color)
set(pump,'color',pump_color)

% energy levels
plot([3.76 6.67],([1 1]*(Eb+wavenum920))*1e-3,'-k','linewidth',1)
plot([3.726 3.989],[1 1]*1*(Eb-wavenum635+wavenum920)*1e-3,'-k','linewidth',1)

hold off;

set(text(5.2,2,'X^1\Sigma^+'),'interpreter','tex','fontsize',10,'color',get(Xplot,'color'));
set(text(7.6,4.1,'a^3\Sigma'),'interpreter','tex','fontsize',10,'color',get(aplot,'color'));
set(text(7.1,15.5,'c^3\Sigma'),'interpreter','tex','fontsize',10,'color',get(cplot,'color'));
set(text(5.1,12.7,'b^3\Pi'),'interpreter','tex','fontsize',10,'color',get(bplot,'color'));
set(text(6.0,16.9,'B^1\Pi'),'interpreter','tex','fontsize',10,'color',get(Bplot,'color'));

set(text(8,5.6,'Na(3s)+Cs(6s)'),'fontsize',10)
set(text(7.8,17.2,'Na(3s)+Cs(6p_{3/2})'),'fontsize',10)
set(text(4.4,0.04965,'v=0, J=0'),'fontsize',10)
set(text(4.,-0.8,'(rovibrational ground state)'),'fontsize',8)

set(text(6.65,8.8,{'\Omega_P','922 nm'}),'fontsize',10,'color',get(pump,'color'))
set(text(4,8.8,{'\Omega_S','635 nm'}),'fontsize',10,'color',get(stokes,'color'))

ylim(([-0.031 0.06]/wavenum2hartree - yoff)/yscale)
xlim([3.6 19]*bohr2angstrom)
set(gca,'fontsize',10,'color','none')

xlabel('Internuclear distance R (A)','interpreter','tex')
ylabel('Energy (10^3 cm^{-1})','interpreter','tex')

set(gcf,'units','inches','position',[0.4861    3.7500    3.75    4]);

end