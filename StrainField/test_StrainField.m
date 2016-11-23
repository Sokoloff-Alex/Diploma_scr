% Test_StrainField 

close all
clear all
clc
%%
clc
x1 = [0 0 0 1 1 1 2 2 2]'+1;
y1 = [0 1 2 0 1 2 0 1 2]'+2;


x2 = x1 + [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.07]';
y2 = y1 + [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.08]'*3;

vel = [x2-x1, y2-y1]/10;


%%
clc
Strain_test = getStrainMap2([x1,y1,vel])

StrainL1 = grid2vector(Strain_test(:,:,1));
StrainL2 = grid2vector(Strain_test(:,:,2));
StrainAngle = grid2vector(Strain_test(:,:,3));
StrainNormal = [StrainL1, StrainL2, StrainAngle];
Strain_Lat = grid2vector(Strain_test(:,:,7));
Strain_Lon = grid2vector(Strain_test(:,:,8));

%%
close all
figure(1)
hold on; grid on; axis equal
plot(x1,y1,'ok')
plot(x2,y2,'or')
quiver(x1,y1, vel(:,1),vel(:,2),'k')
plot(grid2vector(Strain_test(:,:,7)), grid2vector(Strain_test(:,:,8)), '*')
plotStrain(StrainNormal, Strain_Lon, Strain_Lat,10^3);
plot(x1(1:3),y1(1:3),'--k')
plot(x1(4:6),y1(4:6),'--k')
plot(x1(7:9),y1(7:9),'--k')
plot(x2(1:3),y2(1:3),'--r')
plot(x2(4:6),y2(4:6),'--r')
plot(x2(7:9),y2(7:9),'--r')
plot(x1([1,4,7]),y1([1,4,7]),'--k')
plot(x1([2,5,8]),y1([2,5,8]),'--k')
plot(x1([3,6,9]),y1([3,6,9]),'--k')
plot(x2([1,4,7]),y2([1,4,7]),'--r')
plot(x2([2,5,8]),y2([2,5,8]),'--r')
plot(x2([3,6,9]),y2([3,6,9]),'--r')
% xlim([1 5])
% ylim([0 4])