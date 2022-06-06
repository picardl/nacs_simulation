%%%%unit: um

clear all;
clc;

global lamda k n1 NA fo
lamda=623e-3;                               % wavelength
k=2*pi/lamda;                               % propagation (wave) constant free space
n1 = 1;                                     % refractive index
NA = 0.55;                                  % numerical aperture of the objective
Mo = 100;                                   % magnification of the lens
F = 180e3;                                  % tube length
fo=F./Mo;                                   % objective focal length (F=tube Length, M=magnification)
R=fo.*NA./n1;                               % objective back aperture (radius)
d=1000;                                     % propagation distance
w = 10000;                                  % beam waist

% Sampling stuff
L1=2e5; %side length
M=400; %number of samples
dx1=L1/M; %src sample interval
x1=-L1/2:dx1:L1/2-dx1; %src coords
y1=x1;

my0=1080;                                                                   
mx0=1080;                                                        % resolution of the input plane
pixel0=8;                                                        % pixel size
L0=(mx0-1)*pixel0;                                               % dimension of input plane
[xx,yy]=meshgrid(-(my0-1)/2:(my0-1)/2,-(my0-1)/2:(my0-1)/2);
Aperture=sign(1-sign(xx.^2+yy.^2-((my0-1)./2).^2));              % circular aperture
A0=Aperture;
%%%%%% here define the input phase profile 
% g=k/(2*d).*(xx.^2+yy.^2);
g=A0.*exp(-(xx.^2+yy.^2)./w);
% .*exp(1i.*g); % 3-fold helical vortex

figure(1)
imagesc(abs(g.^2));

% Tilting this input beam
deg=pi/180;
alpha=5.0e-5; %rad
theta= 0;%45*deg;

% [u1]=tilt(g,L1,lamda,alpha,theta);
u1 = g;
u2=propFF(u1,L1,lamda,d); %propagation of src 1
I2 = abs(u2.^2);
figure(2) %display obs irrad
imagesc(x1,y1,I2);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title(['Propagated tilted beam; z= ',num2str(d),' m']);