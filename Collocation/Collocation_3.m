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
[CRD, SigmaVenu_merged, name ] = merge_stations(CRD_all,SigmaVenu,names_all);
[CRD, AngleV_merged] = merge_stations(CRD_all,AngleV,names_all)
Angle_v = AngleV_merged(:,1);

%% save Error Bars for GMT
% fileID = fopen('~/Alpen_Check/MAP/VelocityField/Vu_bars.txt', 'w');
% fprintf(fileID, '# Velocity Field Error Bars Lat=Lat+Vu*scale, SigmaVu (mm/yr) -> SigmaVu[deg/yr] (for ploting with "gmt psxy -Ex" ) \n');
% fprintf(fileID, '#  Long [deg],   Lat [deg],      Sigma U [deg/yr],      \n');
% formatStr = '%12.7f  %12.7f   %15e \n';
% d = diag(CovVenu);
% SigmaVu = SigmaVenu_merged(:,3) * 1.8^(1/2) * 20; % scaled to abequate value [m/yr]
% SigmaVu_deg_yr = SigmaVu*1000 * 0.36; 
% data = [long(iiSel), lat(iiSel) + Vu_res(iiSel) * 1000 * 0.24, SigmaVu_deg_yr(iiSel)];
% 
% for i = 1:length(iiSel)
%    fprintf(fileID, formatStr, data(i,:)); 
% end
% fclose(fileID);

%% compute common observation period

t_start = SINEX.SOLUTION.EPOCHS.DATA_START;
t_end   = SINEX.SOLUTION.EPOCHS.DATA_END;

t_start = [str2num(t_start(:,1:2)) + str2num(t_start(:,4:6))/365.25]
t_end   = [str2num(t_end(  :,1:2)) + str2num(t_end(  :,4:6))/365.25]

dt = t_end - t_start;

%%
[dtobs_sum] = merge_stations_sum(dt,names_all);

%%
close all
figure(1)
hist(dtobs_sum,10)

%% save full table of velocity field
% d = diag(CovVenu);
clc
Sigma_Venu = SigmaVenu_merged * 1.8^(1/2) * 20 * 1000;
data = [wrapTo180(long), lat, Ve_res*1000, Sigma_Venu(:,1), Vn_res*1000, Sigma_Venu(:,2), Angle_v, Vu_res*1000, Sigma_Venu(:,3), dtobs_sum];
formatStr = '%4s  %8.3f\t %8.3f\t %5.2f  %5.2f   %5.2f  %5.2f   %5.1f     %5.2f  %5.2f   %4.1f \n';
disp('Site  Long,     Lat,          Ve     SVe     Vn     SVn     angle     Vu     Svu   Tobs')
for i = 1:length(long)
    fprintf(formatStr, name{i}, data(i,:)); 
end


%% for Horizontal 

iiSel = Selected;        
iiOut = setdiff([1:198], iiSel);

%% wrie outliers files
% writeVelocityFieldVerticalOutliers2GMT(long(iiOut), lat(iiOut), Vu_res(iiOut), names(iiOut), '~/Alpen_Check/MAP/VelocityField/VelocityFieldVertical_out.txt');

%%
clc
CovVenu2 = extractCovariance(CovVenu, iiSel, [1 2 3], 'no split');
V_enu_res = [Ve_res(iiSel), Vn_res(iiSel), Vu_res(iiSel)];
%% for Horizontal
clc
%                                                                                                             LongLim LatLing step dMax nMin flags ...                                         
[LongGrid, LatGrid, V_def, rmsFit, V_SigPred] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [-4 17], [41 53], 1, 150, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr');

%% for Vertical       

%% Block selection for Vertical

Outliers   = {'HELM', 'WIEN', 'OGAG', 'OBE2', ...
              'ROHR','BIWI','BI2I', 'MANS'};
            
Outliers = {'ELMO','WIEN','FERR', ...
            'BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE2','WT21','HKBL','PATK','PAT2', ...
            'HRIE','KTZ2', 'WLBH'};
        
iiOut    = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);
[LongGrid, LatGrid, V_def, rmsFit, V_SigPred] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [1 16], [42 49], 1, 250, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%% run Kriging
close all
runKriging(LongGrid, LatGrid, V_def, long ,lat, Vu_res, iiSel);

%% check individual
run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [10 10], [45 45], 1, 250, 10, 'exp1', '-v', 'no bias', 'tail 0', 'no corr', 'filter','plot');

%% Plot Results 2D map
clc
s = 500;  % [mm/yr]
try
    close (fig7)
end
clr = lines(8);
fig7 = figure(7);
hold on
xlim([-2 19])
ylim([41 52])
etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
% Earth_coast(2)
% plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
% plot(Adriatics(:,1),Adriatics(:,2) , '--k')
quiver(long(iiSel), lat(iiSel), Ve_res(iiSel)*s,    Vn_res(iiSel)*s,  0, 'r', 'lineWidth',1)
% quiver(LongGrid,      LatGrid,      V_def(:,1)*s,   V_def(:,2)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% errorbar(long(iiSel), lat(iiSel)+ Vu_res(iiSel)*s, SigmaVenu(iiSel,3)*s, '.m')
% quiver(long(iiSel), lat(iiSel), zeros(size(iiSel))', Vu_res(iiSel)*s,  0, 'r', 'lineWidth',1)
quiver(long(iOutliers),lat(iOutliers),Ve_res(iOutliers)*s,   Vn_res(iOutliers)*s, 0, 'm', 'lineWidth',1)
% quiver(long(iiOut),lat(iiOut),zeros(size(iiOut))',  Vu_res(iiOut)*s, 0, 'm', 'lineWidth',1)
% text(long(iiOut),lat(iiOut), names(iiOut))
quiver(LongGrid,      LatGrid,      V_def(:,1)*s,   V_def(:,2)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% errorbar(LongGrid, LatGrid + V_def(:,3)*s, rmsFit(:,3)*s, '.k')
% errorbar(LongGrid, LatGrid + V_def(:,3)*s, V_SigPred(:,3)*s, '.k')
% quiver(LongGrid,      LatGrid,      zeros(size(LongGrid)),   V_def(:,3)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% text(long,lat, names)
plotErrorElipses('Cov', CovVenu,                         long,     lat,     Ve_res,     Vn_res,     s, 0.95, 'r')
plotErrorElipses('Sig', sqrt([V_SigPred(:,1),V_SigPred(:,2)])*2, LongGrid, LatGrid, V_def(:,1), V_def(:,2), s, 0.95, 'b')
% plotErrorElipses('Sig', ([rmsFit(:,1),rmsFit(:,2)])*1000,  LongGrid, LatGrid, V_def(:,1), V_def(:,2), s, 0.95, 'm')

title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Residual velocity','Outliers','LSC velocity','location','NorthWest')

%%
clc
%  VARIANCE FACTOR                     1.800168208557666
% sc = 20*1.8
% sigmaVenu_SNX = SigmaVenu;
% writeVelocityFieldwithEllipseGMT([ LLVS, AngleV  ], names_all, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor_Ellipse_SNX.txt')
writeDeformationFieldGMT([LongGrid, LatGrid V_def(:,[1,2]), rmsFit(:,[1,2])*1, zeros(size(LongGrid))], '~/Alpen_Check/MAP/VelocityField/VelocityFieldPredictedFitting.txt', '-e')

%%
clc
writeDeformationFieldGMT([LongGrid, LatGrid V_def(:,[1,2]), V_SigPred(:,[1,2]), zeros(size(LongGrid))], '~/Alpen_Check/MAP/VelocityField/VelocityFieldPropagatedNoise.txt', '-e')


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
