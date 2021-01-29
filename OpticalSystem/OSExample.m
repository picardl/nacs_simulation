%% 4f imaging system

OS = OpticalSystem(0, 'mm'); % start z, unit (only a label)
OS.addLens(50,50,25.4); %z,f,diameter (diameter is just for illustration)
OS.addLens(850,750,25.4 *2);
OS.addLens(1600,16,18); % objective

ang = [-1,-1,0,0,1,1] * pi / 180; % technically this is the slope. For small slopes, the slope (tangent of the angle) is approximately the angle
pos = [-1,1,-1,1,-1,1];
colors = {'r','r','b','b','g','g'};

rays = cat(1,pos,ang);

zlist = linspace(0,1600,400);

res = OS.propagateRays(rays, zlist, colors);

%% 2 4f systems
OS = OpticalSystem(0);
OS.addLens(50,50,25.4);
OS.addLens(150,50,25.4);
OS.addLens(250,50,25.4);
OS.addLens(1050,750,25.4 * 2);
OS.addLens(1800,16,18);

ang = [-1,-1,0,0,1,1] * pi / 180; % technically this is the slope. For small slopes, the slope (tangent of the angle) is approximately the angle
pos = [-1,1,-1,1,-1,1];
colors = {'r','r','b','b','g','g'};

rays = cat(1,pos,ang);

zlist = linspace(0,1800,400);

res = OS.propagateRays(rays,zlist,colors);