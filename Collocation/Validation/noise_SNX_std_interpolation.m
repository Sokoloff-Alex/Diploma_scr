% Srcipt to interpolate SNX std values

close all
clear all
clc

%%
SINEX       = readSNX('STA/FMC_IGB_W7.SNX');
ALP_NET_CRD = readCRD('STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('STA/FMC_IGB_W7.VEL');
% [Stations, Radoms, Records] = readOUT('STA/FMC_IGB_W7.OUT');



%%
flags = ALP_NET_CRD(:,7);
range_flag_W = 1:length(flags);
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);
range_flag_A = 1:length(flags);
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);
range_flag = sort([range_flag_A, range_flag_W]);

CRD = cell2mat( ALP_NET_CRD(range_flag,4:6));
VEL = cell2mat( ALP_NET_VEL(range_flag,4:6));

names_all = ALP_NET_CRD(range_flag,2);

Azim_All = Records.VEL.ENU.Ellipse.Angle;
VEL_std = SINEX.SOLUTION.ESTIMATE.Data.VEL_STD;
Sites_list = SINEX.SITE.ID.CODE;
Stations = cellstr(SINEX.SITE.ID.CODE);
DOMES    = cellstr(SINEX.SITE.ID.DOMES);
SiteDome = [SINEX.SITE.ID.CODE, repmat(char(' '),297,1), SINEX.SITE.ID.DOMES];
SiteDome_list = cellstr(SiteDome);

[SigmaVe, SigmaVn, SigmaVu] = XYZ2ENU(CRD, VEL_std);
SigmaVn = -SigmaVn;

% Megre artificial stations
[CRD,VEL,names]   = merge_stations(CRD, VEL,names_all);
[SigmaVenu, Azim] = merge_stations([SigmaVe, SigmaVn, SigmaVu], Azim_All, names_all);
Azim = Azim(:,1);

SigmaVe = SigmaVenu(:,1);
SigmaVn = SigmaVenu(:,2);



% remove Eurasia plate motion
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

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
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
plot(Adriatics(:,1),Adriatics(:,2) , '--k')
% measurements
quiver(long(Selected),lat(Selected),Ve_res(Selected)*sc1,      Vn_res(Selected)*sc1,        0, 'r', 'lineWidth',1)
    sc2 = 50;
    scSigVe = sc1 * sc2 * 100;
    scSigVn = sc1 * sc2 * 350;
for i = 1:length(Selected) 
   ellipce_2D([km2deg(SigmaVe(Selected(i))*scSigVe, 6378*cosd(lat(Selected(i)))), SigmaVn(Selected(i))*scSigVn/111], ...
       Azim(Selected(i)), [long(Selected(i)) + Ve_res(Selected(i))*sc1, lat(Selected(i))  + Vn_res(Selected(i))*sc1], 1, clr(3,:)) 
end
% text(long,lat, names)

% interpolation results
% quiver(LongGrid,      LatGrid,      V_pred_1(:,2)*s,   V_pred_1(:,1)*s,    0, 'Color',clr(2,:), 'lineWidth',1)
quiver(LongGrid,      LatGrid,      V_pred_2(:,2)*sc1,   V_pred_2(:,1)*sc1,    0, 'Color',clr(1,:), 'lineWidth',1)
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


