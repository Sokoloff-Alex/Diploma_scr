% Test_StrainField 

x1 = [0 1 2 0 1 2 0 1 2]';
y1 = [0 0 0 1 1 1 2 2 2]';

x2 = x1 + [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.07]';
y2 = y1 + [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.08]'*3;


x1 = x1 + 1;
x2 = x2 + 1;
y1 = y1 + 1;
y2 = y2 + 1;

vel = [x2-x1, y2-y1];

Strain_test = getStrain([x1,y1,vel]);
Strain_test.NormalStrain
close all
figure(1)
hold on; grid on; axis equal
plot(x1,y1,'ok')
plot(x2,y2,'or')
quiver(x1,y1, vel(:,1),vel(:,2),'k')
plotStrain(Strain_test.NormalStrain,Strain_test.Grid(:,1),Strain_test.Grid(:,2),500*1000);
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
% xlim([-1 2])
% ylim([-1 2])