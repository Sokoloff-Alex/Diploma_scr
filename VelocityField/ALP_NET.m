%% ALP_NET Network
clear all
close all
clc

%%
ALP_NET_CRD = importfileCRDfile('ALP_NET.CRD');
ALP_NET_VEL = importfileVELfile('ALP_NET.VEL');

%%

Marker_CRD = (ALP_NET_CRD(:,2));
X = cell2mat(ALP_NET_CRD(:,4));
Y = cell2mat(ALP_NET_CRD(:,5));
Z = cell2mat(ALP_NET_CRD(:,6));

Vx = cell2mat(ALP_NET_VEL(:,4));
Vy = cell2mat(ALP_NET_VEL(:,5));
Vz = cell2mat(ALP_NET_VEL(:,6));

Rxyz_I = [X,  Y,  Z];
Vxyz_I = [Vx, Vy, Vz];

[Rxyz_E, Vxyz_E] = ITRF2ETRF(Rxyz_I,Vxyz_I);

dv = Vxyz_I + Vxyz_E;

%%

[Long_h,     Lat_h,     velE_h,     velN_h,  Marker] = readVelocityFieldfile('Velocity_field_horizontal.txt');
[Long_Vup,   Lat_Vup,   velE_Vup,   velN_Vup  ]      = readVelocityFieldfile('Velocity_field_vertical_Uplift.txt');
[Long_Vdown, Lat_Vdown, velE_Vdown, velN_Vdown ]     = readVelocityFieldfile('Velocity_field_vertical_subduction.txt');

%%

HELMR1 = readHLMfile('HELMR1.L12',22,230); % NEU
HELMR2 = readHLMfile('HELMR1.L21',22,231); % XYZ

Marker_Res = (HELMR1(:,1));
Res_N = str2num(cell2mat(HELMR1(:,3)));
Res_E = str2num(cell2mat(HELMR1(:,4)));
Res_U = str2num(cell2mat(HELMR1(:,5)));

Marker_Res2 = (HELMR2(:,1));
Res_X = str2num(cell2mat(HELMR2(:,3)));
Res_Y = str2num(cell2mat(HELMR2(:,4)));
Res_Z = str2num(cell2mat(HELMR2(:,5)));

%% 
index1=1;
index2=1;
range=zeros(size(Marker_Res2));
for i = 1:size(Marker_CRD,1)
    if strcmp(Marker_CRD(index1), Marker_Res2(index2))
        range(index2) = i;
        index1=index1 + 1;
        index2=index2 + 1;
    else
        index1=index1 + 1;
    end    
end
clc


%%
clc
close all

figure(1)
hold on
grid on
axis equal
plt1 = plot3(X, Y, Z,'.b');
% plot3(Vx, Vy, Vz,'*b')
plt2 = quiver3(Rxyz_I(:,1),Rxyz_I(:,2),Rxyz_I(:,3),Vxyz_E(:,1),Vxyz_E(:,2),Vxyz_E(:,3),5,'r');
plt3 = quiver3(X(range),Y(range),Z(range),Res_X,Res_Y,Res_Z,10);
%quiver3(Rxyz_I(:,1),Rxyz_I(:,2),Rxyz_I(:,3),Vxyz_I(:,1),Vxyz_I(:,2),Vxyz_I(:,3),'r')
% quiver3(Rxyz_I(:,1),Rxyz_I(:,2),Rxyz_I(:,3),dv(:,1),dv(:,2),dv(:,3),'g')
plt4 = text(X, Y, Z,Marker_CRD);
% Earth_coast(3)
% xlim([min(X) max(X)])
% ylim([min(Y) max(Y)])
% zlim([min(Z) max(Z)])
xlabel('X ,[m]')
ylabel('Y ,[m]')
zlabel('Z ,[m]')
legend([plt1, plt2, plt3]','Station', 'Velocity', 'Residuals')



%%
close all
fig1 = figure(1);
hold on
title('Velocity field map');
xlabel('Longitude, [deg]');
ylabel('Latitude, [deg]');
vhor   = quiver(Long_h,     Lat_h,     velN_h,     velE_h , 'b');
vvdown = quiver(Long_Vdown, Lat_Vdown, velN_Vdown, velE_Vdown, 'r');
vvup   = quiver(Long_Vup,   Lat_Vup,   velN_Vup,   velE_Vup, 'g');
legend([vhor vvdown vvup]','Horisontal velocity','Vertical velocity Up', 'Vertical velocity down');
text(Long_h, Lat_h, Marker);
% vhor = quiver(Long_h, Lat_h, Res_E, Res_N , 'm');
% Earth_coast(2)
xlim([-10 20])
ylim([40 55])
set
hold off

%%

