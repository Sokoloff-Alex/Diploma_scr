% Srcipt to interpolate SNX std values

close all
clear all
clc

%%
[Stations, Radoms, Records] = readOUT('STA/FMC_IGB_W7.OUT');

%%
SINEX       = readSNX('STA/FMC_IGB_W7.SNX','All');
ALP_NET_CRD = readCRD('STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('STA/FMC_IGB_W7.VEL');


%%
flags = ALP_NET_CRD(:,7);
range_flag_W = 1:length(flags);
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);
range_flag_A = 1:length(flags);
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);
range_flag = sort([range_flag_A, range_flag_W]);

CRD_all = cell2mat( ALP_NET_CRD(range_flag,4:6));
VEL_all = cell2mat( ALP_NET_VEL(range_flag,4:6));

names_all = ALP_NET_CRD(range_flag,2);
DOMES = ALP_NET_CRD(range_flag,3);

% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD_all,VEL_all,names_all);

% remove Eurasia plate motion
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);
[V_res_xyz_all] = remove_plate_motion(CRD_all, VEL_all, Omega_Eur);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]
[Ve_res_all, Vn_res_all, Vu_res_all] = XYZ2ENU(CRD_all,V_res_xyz_all); % NEU components, [m/yr m/yr m/yr]

[VEL_std_xyz, VEL_std_xyz2] = merge_stations(VEL_std, VEL_std, names_all);

%%
Azim_All = Records.VEL.ENU.Ellipse.Angle;
VEL_std = SINEX.SOLUTION.ESTIMATE.Data.VEL_STD;

Sites_list = SINEX.SITE.ID.CODE;
Stations = cellstr(SINEX.SITE.ID.CODE);
DOMES    = cellstr(SINEX.SITE.ID.DOMES);
SiteDome = [SINEX.SITE.ID.CODE, repmat(char(' '),297,1), SINEX.SITE.ID.DOMES];
SiteDome_list = cellstr(SiteDome);

%% Transform sigma / covariance from XYZ to ENU

clc
SNX_cov = SINEX.SOLUTION.COVA_ESTIMATE;

close all
figure(2)
hold on; grid on
axis equal
sigmaVenu = zeros(size(VEL_std_xyz));
angle = size(VEL_std_xyz,1);
corrVen = size(VEL_std_xyz,1);
for i = 1:297 %198 %size(VEL_std_xyz,1)
%     covVxyz = [VEL_std_xyz(i,1)^2 0 0 ;
%                0 VEL_std_xyz(i,2)^2 0 ;
%                0 0 VEL_std_xyz(i,3)^2];
%    [covVenu, sigmaVenu(i,:)] = covXYZ2ENU(covVxyz, lat(i), long(i));

    iCv = [((i-1)*6+4) : (((i-1)*6+6))];
    covVxyzSNX = SNX_cov(iCv,iCv);
    covVxyzSNX(1,2) = covVxyzSNX(2,1);
    covVxyzSNX(1,3) = covVxyzSNX(3,1);
    covVxyzSNX(2,3) = covVxyzSNX(3,2);
    [covVenu, sigmaVenu(i,:)] = covXYZ2ENU(covVxyzSNX, lat_all(i), long_all(i));
    %     SigmaVxyz_Cov(i,1:3) = sqrt( [covVxyzSNX(1,1), covVxyzSNX(2,2), covVxyzSNX(3,3)] );
    SigmaVenu_Cov(i,1:3) = sqrt( [covVenu(1,1), covVenu(2,2), covVenu(3,3)] );
    corrVen(i,:) = covVenu(1,2) / (sigmaVenu(i,1) * sigmaVenu(i,2));
    angle(i,:) = 90 + 1/2 * atand(2*covVenu(1,2) / (covVenu(1,1) - covVenu(2,2)));
    if sigmaVenu(1) < sigmaVenu(1);
        angle(i,:) = angle(i,:) + 90;
    end
    error_ellipse([covVenu(1:2, 1:2)]*1000^2, [0 0 ],0.67, 10, 'r')
end
azim = 90 - angle;

%%
%  VARIANCE FACTOR                     1.800168208557666
sc = 10*1.8
% sigmaVenu_SNX = sigmaVenu;
LLVS = [long_all , lat_all, Ve_res_all, Vn_res_all, sigmaVenu_SNX(:,1)*1.8*30, sigmaVenu_SNX(:,2)*1.8*20];
writeVelocityFieldwithEllipseGMT([ LLVS, angle  ], names_all, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Ellipse_SNX.txt')
writeVelocityFieldwithCovGMT(    [ LLVS, corrVen], names_all, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Cov_SNX.txt')

%%
sc = 10
% sigmaVenu_std = sigmaVenu;
LLVS = [long , lat, Ve_res, Vn_res, sigmaVenu_std(:,1)*sc*1, sigmaVenu_std(:,2)*sc];
writeVelocityFieldwithEllipseGMT([ LLVS, angle  ], names, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Ellipse_std.txt')
writeVelocityFieldwithCovGMT(    [ LLVS, corrVen], names, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Cov_std.txt')


fileID = fopen('~/Alpen_Check/MAP/VelocityField/list_outliers','w');
fprintf(fileID, '%4s\n', Outliers{:});
fclose(fileID);



%%
clear V_pred_2 Sigma_V_pred_2 LatGrid LongGrid
tic
range = 1:length(lat);
Max_Dist = 150; % km
lim = 10;
p = 0;
step = 1;
for iLong = -4:step:18
    for iLat = 42:step:53
        arc = greatcircleArc(iLat, iLong, lat, long) * 111 ; % km
        sel = range(arc < Max_Dist);    
        sel = intersect(sel, Selected);
        if length(sel) < lim % add more stations
            add = sort(arc);
            add = add(1:lim); % ad 2 sites, ast is prop already included
            iadd = ismember(arc,add);
            inew = range(iadd);
            sel = unique(sort([sel, inew]));
            sel = intersect(sel, Selected);
        end
        if length(sel) > 1 % do LSC !!! 
            p = p + 1; % point Number
%             V_pred_1(p,:)       = solve_LSC(iLat, iLong, lat(sel),
%             long(sel), Vn_res(sel)*1000,  Ve_res(sel)*1000, 'exp1', '-v', 'bias', 'tail 0')'/1000;   % LSC, no diff, exept few outliers 
%             V_pred_2(p,:)       = solve_LSC(iLat, iLong, lat(sel), long(sel), Vn_res(sel)*1000,  Ve_res(sel)*1000, 'exp1', '-v', 'bias', 'tail 0', 'no corr')'/1000;   % LSC
            V_pred_2(p,:) = solve_LSC_2(iLat, iLong, lat(sel), long(sel), Vn_res(sel)*1000,  Ve_res(sel)*1000, SigmaVn(sel)*1000, SigmaVe(sel)*1000,'exp1','-v', 'bias', 'tail 0', 'no corr')'/1000;   % LSC *1000/1000 to avoid singularity  
            LatGrid(p,1)  = iLat;
            LongGrid(p,1) = iLong;
        end
    end
end
t2 = toc

% correlation error
CorrErr = V_pred_1 - V_pred_2;
Sigma_V_pred_2 = V_pred_2(:,[3,4]);

%% Plot Results 2D map
clc
sc1 = 500;  % [mm/yr]
try
    close (fig7)
end
clr = lines(8);
fig7 = figure(7);
hold on

% xlim([2 16])
% ylim([43 49])
% geoshow(Etopo_Europe, refvec_Etopo, 'DisplayType', 'texturemap');
% demcmap(Etopo_Europe);
% cptcmap('Europe')
% Earth_coast(2)
plot(wrapTo180(Orogen_Alp(:,1)),Orogen_Alp(:,2),'--m')
plot(wrapTo180(Adriatics(:,1)),Adriatics(:,2) , '--k')
% measurements
quiver(wrapTo180(long(Selected)),lat(Selected),Ve_res(Selected)*sc1,      Vn_res(Selected)*sc1,        0, 'r', 'lineWidth',1)
    sc2 = 1;
    scSigVe = sc1 * sc2 * 100;
    scSigVn = sc1 * sc2 * 350;
for i = 1:length(Selected) 
   ellipce_2D([km2deg(SigmaVe(Selected(i))*scSigVe, 6378*cosd(lat(Selected(i)))), SigmaVn(Selected(i))*scSigVn/111], ...
       Azim(Selected(i)), [long(Selected(i)) + Ve_res(Selected(i))*sc1, lat(Selected(i))  + Vn_res(Selected(i))*sc1], 1, clr(3,:)) 
end
% text(long,lat, names)

% interpolation results
% quiver(LongGrid,      LatGrid,      V_pred_1(:,2)*s,   V_pred_1(:,1)*s,    0, 'Color',clr(2,:), 'lineWidth',1)
quiver(LongGrid,      LatGrid,      V_def(:,2)*sc1,   V_def(:,1)*sc1,    0, 'Color',clr(1,:), 'lineWidth',1)
% quiver(LongGrid,      LatGrid,        CorrErr(:,2)*s,    CorrErr(:,1)*s,     0, 'Color',clr(4,:), 'lineWidth',1)
for i = 1:length(LongGrid) 
   ellipce_2D([km2deg(Sigma_V_pred_2(i,2)*scSigmaVe,6378*cosd(LatGrid(i))), Sigma_V_pred_2(i,1)*scSigmaVn/111], ...
0, [LongGrid(i) + V_pred_2(i,2)*s, LatGrid(i) + V_pred_2(i,1)*s], 1) 
end

% ellipce_2D([km2deg(Max_Dist,6378*cosd(iLat)), Max_Dist/111], 0, [iLong, iLat], 1)
title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Earth Coast','Alps Orogen boundary','Ardiatics','Residual velocity','LSC velocity','location','NorthWest')
% impoly(gca, [Adriatics])
xlim([-6 19])
ylim([41 53])

hold off


