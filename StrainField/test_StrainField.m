% Test_StrainField 

close all
clear all
clc
%% test 1
clc
x =  [0 0 0 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12 13 13 13 14 14 14 15 15 15 16 16 16 17 17 17 18 18 18]';
y =  [4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6 4 5 6  4  5  6  4  5  6  4  5  6  4  5  6  4  5  6  4  5  6  4  5  6  4  5  6  4  5  6]';
vx = [0 0.01 0.01 0.01 0.02 0.03 0.05 0.05 0.07 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.1 0.1 -.1 -.1 -.1 0.1 0.1 0.1 0.1 0.1 0.1 -.1 -.1 -.1 0.0 0.0 0.0 0.0 0.0 0.0 0.1 0.1 0.1 -.1 -.1 -.1 0.1 0.1 0.1]'*1/1000;
vy = [0 0.01 0.01 0.01 0.02 0.04 0.05 0.05 0.08 0.1 0.2 0.5 0.1 0.2 0.3 0.3 0.2 0.2 0.3 0.2 0.0 0.3 0.2 0.0 0.6 0.2 0.0 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 -.1 -.1 -.1 -.1 -.1 -.1 -.0 -.1 -.2 -.0 -.1 -.2 0.1 -.1 0.1 0.1 -.1 0.1 0.1 -.1 0.1]'*1/1000;

Deformation = [x, y, vx, vy];
%% test 2
% clc
% n = 40;
% x = repmat([1:n], n,1);
% y = repmat([1:n]',1,n);
% 
% n = 20;
% vx = [repmat(sind([-180:45:180]), 9,1), ones(9,n-9-5), -ones(9,5); ones(n-9,n-5), -ones(n-9,5)  ] / 1000;
% vy = [repmat(sind([-180 -135 -90 -45 0 45 45 90 90 135 135 180 180 45 45 -45 -45 0 0 90]), n,1)']  / 1000;
% 
% vx = [vx, -vx; 
%       vx, -vx];
% vy = [vy, -vy;
%      -vy,  vy];
% 
% Deformation = [grid2stack(x),grid2stack(y),grid2stack(vx1),grid2stack(vy1)];

%% test 3 
clc
n = 8;
x = repmat([1:n], n,1);
y = repmat([1:n]',1,n);

vx = [0  0 -1  1 -1  1  1  1
      0  0  1  1 -1  1 -1  1
      0  0  1  1 -1  1  1 -1
      0  0  1  1 -1  1  1  0
      0  0 -1  1 -1  1  1  0
      0  0  1  1  0 -1  1 -1
      0  0  1  1  0 -1  1 -1
      0  0  0  0  0  0  0  0] /(1000);
  
vy = [0  0 -1  1  1  1 -1 -1
      1  1  1  1  1  1  1  1
      0  0 -1 -1 -1 -1  1  0
     -1 -1  1  1  1  1  0  0
     -1 -1  1  1  1  1  0  0
      1  1  0  0  0  0  0  0
      0  0  0  0  0  0  0  0
      0  0  0  0  0  0  0  0] / (1000);

Deformation = [grid2stack(x), grid2stack(y), grid2stack(vx), grid2stack(vy)];
  

%%

StrainStack= getStrainMap2(Deformation)

%%

writeStrain2GMT(StrainStack, '~/Alpen_Check/MAP/Strain/StrainTEST.txt')

%%
sc = 300;
close all
figure(1)
hold on; grid on; 
axis equal
plot(x,y,'ok')
% plot(x+vx,y+vy,'or')
quiver(grid2stack(x),grid2stack(y), grid2stack(vx)*sc,grid2stack(vy)*sc,0,'k')
plotStrainNormal(StrainStack,10^7);
% xlim([15 20])
% ylim([3 7])

