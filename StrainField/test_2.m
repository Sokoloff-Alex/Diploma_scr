clear all
close all
clc

%%
clc

x1 = [0   0.0     750.004];
y1 = [0 999.915    -0.016];

x2 = [0   0.124  749.994 ];
y2 = [0 999.998   -0.047 ];

vel = [x2-x1; y2-y1]';

pointA = [x1(1) y1(1)];
pointB = [x1(2) y1(2)];
pointC = [x1(3) y1(3)];

velA = vel(1,:);
velB = vel(2,:);
velC = vel(3,:);

Strain_test2 = getStrain4(pointA, pointB, pointC, velA, velB, velC)

%%
close all
figure(1)
hold on
grid on
axis equal
plot(x1, y1, 'ok')
s = 1000;
quiver(x1', y1', vel(:,1)*s, vel(:,2)*s,0, 'm')
plotStrainNormal(Strain_test2,10^6)
text(x1-50, y1+35, {'A', 'B','C'})


