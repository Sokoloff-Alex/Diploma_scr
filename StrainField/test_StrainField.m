% Test_StrainField 

close all
clear all
clc
%%
clc
x = [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11]';
y = [4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6  4  5  6  4  5  6]';
vx = [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.07 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.2 0.2 0.2 -.2 -.2 -.2 0.2 0.2 0.2]'*1/1000;
vy = [0 0.01 0.01 0.01 0.02 0.04 0.05 0.05 0.08 0.1 0.2 0.5 0.1 0.2 0.5 0.2 0.2 0.2 0.3 0.2 0.0 0.3 0.2 0.0 0.6 0.2 0.0 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2 0.2]'*2/1000;

vel = [vx, vy];
Deformation = [x,y,vel];

StrainStack= getStrainMap2(Deformation)

%%
close all
figure(1)
hold on; grid on; axis equal
plot(x,y,'ok')
% plot(x+vx,y+vy,'or')
quiver(x,y, vel(:,1),vel(:,2),'k')
% plot(grid2vector(Strain_test(:,:,7)), grid2vector(Strain_test(:,:,8)), '*')
% plot(StrainStack(:,1),StrainStack(:,2), 'x')
plotStrain(StrainStack,10^8);

