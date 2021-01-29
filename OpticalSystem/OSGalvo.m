%% 2 4f systems
OS = OpticalSystem(0, 'mm'); %z,f,diameter (diameter is just for illustration)
OS.addLens(80,80,25.4);
OS.addLens(205,125,25.4);
OS.addLens(305,100,25.4);
OS.addLens(905,500,25.4 * 2);
OS.addLens(1405,18,25.4); % objective

ang = (1/6)*[-15e-6,-15e-6,15e-6,15e-6]; % technically this is the slope. For small slopes, the slope (tangent of the angle) is approximately the angle
pos = [-1,1,-1,1];
colors = {'r','r','b','b','g','g'};

rays = cat(1,pos,ang);

zlist = linspace(0,1430,400);

res = OS.propagateRays(rays,zlist,colors);
ylim([-2e-4,2e-4])
xlim([1423.184,1423.188])