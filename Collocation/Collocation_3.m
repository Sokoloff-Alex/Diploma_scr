% Collocation
% 
% make collocation procedure to make deformation field from velocity field
% do in 3D, use Vu_pred furter in kriging 
% Alexandr Sokolov, KEG
% 07.12.2016

%%
clear all
close all
clc

%% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Orogen_Alp  = importOrogen('dat/PB2002_orogen_Alps.txt');
lwmask = struct2array(load('dat/lwmask25.mat'));

% [Etopo_Europe, refvec_Etopo] = etopo('../../../../MAP/etopo1_bed_c_f4/etopo1_bed_c_f4.flt', 1, [40 54], [-7 19]); % ETOPO
load('ETOPO_Alps.mat')
Adriatics = struct2array(load('dat/Adriatics.mat'));
%
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
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

%% Covariance

load('dat/SINEX.mat')
SNX_cov = SINEX.SOLUTION.COVA_ESTIMATE;

[CovVenuSNX, SigmaVenu, CorrVen, AngleV] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'VEL');
[CovRenuSNX, SigmaRenu, CorrRen, AngleR] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'CRD');

[CovVenu] = megreCov(CovVenuSNX, names_all);

%% Block selection

Outliers   = {'HELM', 'WIEN', 'OGAG', 'OBE2', ...
              'ROHR','BIWI','BI2I', 'MANS'};
            
Outliers = {'HELM','WIEN','FERR', ...
            'BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE2','WT21','HKBL','PATK','PAT2', ...
            'HRIE','KTZ2', 'WLBH'};
   
iiOut    = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);

%% Estimate and remove block motion from ITRF
% [Vel_res_Blocks, Vel_blocks, Omega_Blocks]= remove_block_motion(CRD, VEL, Sets);
% 
% %% Remome Plate motion from Block motion
% [V_bl_no_plate] = remove_plate_motion(CRD, Vel_blocks, Omega_Eur);
% [Ve_res_bl, Vn_res_bl, Vu_res_bl] = XYZ2ENU(CRD,V_bl_no_plate); 
% dVe = Ve_res - Ve_res_bl;
% dVn = Vn_res - Vn_res_bl;
% t = toc;
% disp(['proc time : ', num2str(t), ' [sec]'])

%% LSC, % allocation domain for grid points
clear LatGrid LongGrid V_pred V_pred* V_def*
clc
tic
range = 1:length(lat);
Max_Dist = 250; % km
lim = 10;
p = 0;
step = 0.5;
for iLong = -4:step:18
    for iLat =  42:step:53
        arc = greatcircleArc(iLat, iLong, lat, long) * 111 ; % km
        sel = range(arc < Max_Dist);    
        sel = intersect(sel, iiSel);
        if length(sel) < lim % add more stations
            add = sort(arc);
            add = add(1:lim); % ad 2 sites, ast is prop already included
            iadd = ismember(arc,add);
            inew = range(iadd);
            sel = unique(sort([sel, inew]));
            sel = intersect(sel, iiSel);
        end
        if length(sel) > 1 % do LSC !!! 
            p = p + 1; % point Number
            CovVenuSel = extractCovariance(CovVenu, sel, [1 2 3], 'split');
            Venu = [Ve_res(sel), Vn_res(sel), Vu_res(sel)]*1000;
            [V_pred, rmsFitting, V_noise_pred] = solve_WLSC3(iLat, iLong, lat(sel), long(sel), Venu ,CovVenuSel,'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');   % WLSC 
            V_def3(p,:) = V_pred/1000;
            rmsFit(p,:)   = rmsFitting/1000; % [mm/yr]
            V_SigPred(p,:) = V_noise_pred;
            LatGrid(p,1)  = iLat;
            LongGrid(p,1) = iLong;
        end
    end
end
t2 = toc


%% Plot Results 2D map
clc
s = 250;  % [mm/yr]
try
    close (fig7)
end
clr = lines(8);
fig7 = figure(7);
hold on
xlim([-6 19])
ylim([41 53])
% xlim([2 16])
% ylim([43 49])
% geoshow(Etopo_Europe, refvec_Etopo, 'DisplayType', 'texturemap');
% demcmap(Etopo_Europe);
% cptcmap('Europe')
etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
% Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
plot(Adriatics(:,1),Adriatics(:,2) , '--k')
% quiver(long(Selected), lat(Selected), Ve_res(Selected)*s,    Vn_res(Selected)*s,  0, 'r', 'lineWidth',1)
quiver(long(iiSel), lat(iiSel), zeros(size(iiSel))', Vu_res(iiSel)*s,  0, 'r', 'lineWidth',1)
% quiver(long(iOutliers),lat(iOutliers),Ve_res(iOutliers)*s,   Vn_res(iOutliers)*s, 0, 'm', 'lineWidth',1)
quiver(long(iiOut),lat(iiOut),zeros(size(iiOut))',  Vu_res(iiOut)*s, 0, 'm', 'lineWidth',1)
text(long(iiOut),lat(iiOut), names(iiOut))
% quiver(LongGrid,      LatGrid,      V_def3(:,1)*s,   V_def3(:,2)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
quiver(LongGrid,      LatGrid,      zeros(size(LongGrid)),   V_def3(:,3)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% ellipce_2D([km2deg(Max_Dist,6378*cosd(iLat)), Max_Dist/111], 0, [iLong, iLat], 1)
% text(long,lat, names)
title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Earth Coast','Alps Orogen boundary','Ardiatics','Residual velocity','LSC velocity','location','NorthWest')

%% 

sc1 = 1000;
scSigVe = sc1;
scSigVn = sc1;
for i = 1:length(iiSel) 
    az = Azim(iiSel(i));
    sigma_ee = SigmaVe(iiSel(i)) * scSigVe;
    sigma_nn = SigmaVn(iiSel(i)) * scSigVn;
    sigma_en = -1/2 * tand(2*az) * abs((sigma_ee^2 - sigma_nn^2));
    corrVen(i,1) = sigma_en / (sigma_ee * sigma_nn);
    cov_mat = [sigma_ee^2, sigma_en   ; 
               sigma_en,   sigma_nn^2];
    mean_value = [long(iiSel(i)) + Ve_res(iiSel(i))*sc1, lat(iiSel(i))  + Vn_res(iiSel(i))*sc1];
%    ellipce_2D([sigma_xx sigma_yy], 0,mean_value , 1, clr(3,:)) 
%    error_ellipse(cov_mat, mean_value, 0.95, 'r') % 2 sigma, 95 % confidence
end
grid on
hold off

% Alps_deformation = [LongGrid, LatGrid, V_def1(:,2), V_def1(:,1)];

%%
% clc
% for i = 1:length(SigmaVe) 
%     az = Azim(i);
%     sigma_ee = SigmaVe(i);
%     sigma_nn = SigmaVn(i);
%     sigma_en = 1/2 * tand(2*(90-az)) * abs((sigma_ee^2 - sigma_nn^2));
%     corrVen(i,1) = sigma_en / (sigma_ee * sigma_nn);
% end
% 
% [SigmaVe, SigmaVn, Azim, corrVen]
%%
% sc = 10
% writeVelocityFieldwithEllipseGMT([long(Selected) , lat(Selected), Ve_res(Selected), Vn_res(Selected), SigmaVe(Selected)*100/2, SigmaVn(Selected)*100/2, Azim(Selected)], names(Selected), '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Ellipse.txt')
% writeVelocityFieldwithCovGMT(    [long(Selected) , lat(Selected), Ve_res(Selected), Vn_res(Selected), SigmaVe(Selected)*100/2, SigmaVn(Selected)*100/2, corrVen(Selected)], names(Selected), '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Cov.txt')

% %% save results
% Alps_deformation = [LongGrid, LatGrid, V_def1(:,2), V_def1(:,1)];
% name = 'Alps_deformation_0.25x0.25_no_correlation_3';
% save([name, '.mat'],'Alps_deformation');
% 
% fileID = fopen([name,'.txt'], 'w');
% headString = '  Long [deg],   Lat [deg],  Vel E [m/yr],  Vel N [m/yr] \n';
% formatStr = '%12.7f  %12.7f  %12.5f  %12.5f \n';
% fprintf(fileID, 'Deformation / Velocity field horizontal \n');
% fprintf(fileID, headString);
% for i = 1:length(LongGrid)
%    fprintf(fileID, formatStr, Alps_deformation(i,:)); 
% end
% fclose(fileID);
