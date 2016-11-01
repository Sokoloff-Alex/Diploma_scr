% ALP_EPN_A solution

clear all
close all
clc

%% load data
CRD = importfileCRDfile('../Results/FCSEPN.CRD');
VEL = importfileVELfile('../Results/FCSEPN.VEL');

CRD = importfileCRDfile('../Results/FCSIGSB.CRD');
VEL = importfileVELfile('../Results/FCSIGSB.VEL');
%%
CRD = importfileCRDfile('../Results/FCS_outl.CRD');
VEL = importfileVELfile('../Results/FCS_outl.VEL');

% CRD = importfileCRDfile('../Results/Fd2CSIGSA.CRD');
% VEL = importfileVELfile('../Results/Fd2CSIGSA.VEL');


%% flag filter 
flags = (CRD(:,7));

range_flag_A = [1:length(flags)];
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);

range_flag_W = [1:length(flags)];
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);

range_flag_All = [range_flag_A, range_flag_W];
range_flag_All = sort(range_flag_All);



clc
Marker_CRD = (CRD(:,2));
X = cell2mat(CRD(:,4));
Y = cell2mat(CRD(:,5));
Z = cell2mat(CRD(:,6));

Vx = cell2mat(VEL(:,4));
Vy = cell2mat(VEL(:,5));
Vz = cell2mat(VEL(:,6));

R_I = [X,  Y,  Z];
V_I = [Vx, Vy, Vz];

[R_E, V_E]   = ITRF2ETRF(R_I,V_I);
%%
[Ve, Vn, Vu, lat, lon, h] = XYZ2ENU([X, Y, Z] ,[Vx, Vy, Vz]);
[V_vertical] = VerticalProjection(X, Y, Z ,Vx, Vy, Vz);

counterS = 0;
counterU = 0;
for i=1:size(V_vertical)
   V_mag(i) = norm(V_vertical(i,:));
   e_r   = R_I(i,1:3)./norm(R_I(i,1:3));
   e_Vup = V_vertical(i,:)./norm(V_vertical(i,:));
   if norm(e_r + e_Vup) > 1      && strcmp(flags(i), 'A') || strcmp(flags(i), 'W')  && V_mag(i) < 0.004    
       counterU = counterU + 1 ;
       Flag_U(counterU) = i;
   elseif norm(e_r + e_Vup) < 1  && strcmp(flags(i), 'A') || strcmp(flags(i), 'W')  && V_mag(i) < 0.004
       counterS = counterS + 1 ;
       Flag_S(counterS) = i; 
   end
end

Flags_All_good = sort([Flag_U, Flag_S]);

%%
close all
clc
clr = lines(5);
scale = 1000*1000*100;
figure(1)
hold on
grid on
axis equal
% plot3(0,0,0,'*r')
plt1 = plot3(R_E(range_flag_All,1),R_E(range_flag_All,2),R_E(range_flag_All,3),'.b');
% plt2 = quiver3(R_E(:,1),R_E(:,2),R_E(:,3),V_E(:,1),V_E(:,2),V_E(:,3),10,'r');
% plt3 = quiver3(R_E(range_flag_W,1),  R_E(range_flag_W,2),  R_E(range_flag_W,3),  V_E(range_flag_W,1)*scale,  V_E(range_flag_W,2)*scale,  V_E(range_flag_W,3)*scale,  0,'b');
% plt4 = quiver3(R_E(range_flag_A,1),  R_E(range_flag_A,2),  R_E(range_flag_A,3),  V_E(range_flag_A,1)*scale,  V_E(range_flag_A,2)*scale,  V_E(range_flag_A,3)*scale,  0,'r');
% plt5 = quiver3(R_E(range_flag_All,1),R_E(range_flag_All,2),R_E(range_flag_All,3),V_E(range_flag_All,1)*scale,V_E(range_flag_All,2)*scale,V_E(range_flag_All,3)*scale,0,'r');
% 
% plt6 = quiver3(R_I(:,1),R_I(:,2),R_I(:,3),V_I(:,1)*scale, V_I(:,2)*scale, V_I(:,3)*scale,0,'m')
% plt6 = quiver3(R_I(range_flag_All,1),R_I(range_flag_All,2),R_I(range_flag_All,3),V_Vertical(range_flag_All,1)*scale, V_Vertical(range_flag_All,2)*scale, V_Vertical(range_flag_All,3)*scale,0,'Color',clr(1,:),'LineWidth', 2)
plt6 = quiver3(R_I(range_flag_W,1),R_I(range_flag_W,2),R_I(range_flag_W,3),V_vertical(range_flag_W,1)*scale, V_vertical(range_flag_W,2)*scale, V_vertical(range_flag_W,3)*scale,0,'Color',clr(1,:),'LineWidth', 2)
plt7 = quiver3(R_I(Flag_U,1),R_I(Flag_U,2),R_I(Flag_U,3),V_vertical(Flag_U,1)*scale, V_vertical(Flag_U,2)*scale, V_vertical(Flag_U,3)*scale,0,'Color',clr(2,:),'LineWidth', 2)
plt8 = quiver3(R_I(Flag_S,1),R_I(Flag_S,2),R_I(Flag_S,3),V_vertical(Flag_S,1)*scale, V_vertical(Flag_S,2)*scale, V_vertical(Flag_S,3)*scale,0,'Color',clr(3,:),'LineWidth', 2)


% plt6 = quiver3(R_I(:,1),R_I(:,2),R_I(:,3),V_hor(:,1)*scale, V_hor(:,2)*scale, V_hor(:,3)*scale,0,'b')
%  plt8 = text(R_E(range_flag_All,1),R_E(range_flag_All,2),R_E(range_flag_All,3), Marker_CRD(range_flag_All));
% Earth_coast(3)
% xlim([min(X) max(X)])
% ylim([min(Y) max(Y)])
% zlim([min(Z) max(Z)])
title('Velocity field')
xlabel('X ,[m]')
ylabel('Y ,[m]')
zlabel('Z ,[m]')
legend('Station', 'Velocity W', 'Velocity A, Uplift','Velocity A, Subduction')

%%
close all


[xq,yq] = meshgrid(-5:0.5:20, 40:0.5:52);
Flags_All_good = Flags_All_good(abs(Vu(Flags_All_good)) <0.003);
Vu(Vu(:) > 0.003) = 0.003;
Vu(Vu(:) < -0.003) = -0.003;

s = 1000;

figure
% subplot(1,2,1)
vq = griddata(lon(Flags_All_good),lat(Flags_All_good),Vu(Flags_All_good)*1000,xq,yq, 'natural');
meshc(xq,yq,vq);
hold on
plot3(lon(Flags_All_good),lat(Flags_All_good),Vu(Flags_All_good)*1000,'o');
quiver3(lon(Flag_U),lat(Flag_U),zeros(size(Vu(Flag_U))),zeros(size(Vu(Flag_U))),zeros(size(Vu(Flag_U))), Vu(Flag_U)*s,0)
quiver3(lon(Flag_S),lat(Flag_S),zeros(size(Vu(Flag_S))),zeros(size(Vu(Flag_S))),zeros(size(Vu(Flag_S))), Vu(Flag_S)*s,0)
quiver3(lon(range_flag_W),lat(range_flag_W),zeros(size(Vu(range_flag_W))),zeros(size(Vu(range_flag_W))),zeros(size(Vu(range_flag_W))), Vu(range_flag_W)*s,0)
% hcb = colorbar;
% ylabel(hcb,'[mm/yr]')
% caxis([-3 3])
title('Interpolated velocities (estimated and reference)')
xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
zlabel('Vertical velocity [mm/yr]')
xlim([-5 20])
ylim([40 54])
% 
% subplot(1,2,2)
% vq2 = griddata(lon(range_flag_W),lat(range_flag_W),Vu(range_flag_W)*1000,xq,yq, 'natural');
% meshc(xq,yq,vq2);
% hold on
% plot3(lon(range_flag_W),lat(range_flag_W),Vu(range_flag_W)*1000,'o');
% quiver3(lon(Flag_U),lat(Flag_U),zeros(size(Vu(Flag_U))),zeros(size(Vu(Flag_U))),zeros(size(Vu(Flag_U))), Vu(Flag_U)*s,0)
% quiver3(lon(Flag_S),lat(Flag_S),zeros(size(Vu(Flag_S))),zeros(size(Vu(Flag_S))),zeros(size(Vu(Flag_S))), Vu(Flag_S)*s,0)
% quiver3(lon(range_flag_W),lat(range_flag_W),zeros(size(Vu(range_flag_W))),zeros(size(Vu(range_flag_W))),zeros(size(Vu(range_flag_W))), Vu(range_flag_W)*s,0)
hcb = colorbar;
ylabel(hcb,'[mm/yr]')
caxis([-3 3])
% title('Interpolated velocities (reference only)')
% xlabel('Longitude, [deg]')
% ylabel('Latitude, [deg]')
% zlabel('Vertical velocity [mm/yr]')
% xlim([-5 20])
% ylim([40 54])

%%
SigmaFile=readSigmaFile('Results/SSC_CRD.SIG');
CovMat = zeros(size(SigmaFile,1),3,3);
%%
close all
figure(1)
axis equal
hold on
for i = 1:size(SigmaFile,1)
    SigmaX = cell2mat((SigmaFile(i,3)));
    SigmaY = cell2mat((SigmaFile(i,4)));
    SigmaZ = cell2mat((SigmaFile(i,5)));
    
    CovMat_Stack(i,:,:) = [SigmaX   0     0; 
                             0     SigmaY  0;
                             0       0   SigmaZ]^2;
    covMat = squeeze(CovMat_Stack(i,:,:)) + 1e-6*eye(3);     
    error_ellipse(covMat*100, [rand(3,1) ])
                  
end

%%

clc

parseSNX('../Results/Solution_Estimate_Table')

