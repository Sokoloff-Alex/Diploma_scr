% combined multiyear ALP_NET + SAPOS = ALP_SAP
% from weekly solutions

%% SAPOS network
%
% 
clc
close all
clear all
%% read Solution files
SNX_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/SOL/FMC_C_W1.SNX'
CRD_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/STA/FMC_C_W1.CRD'
VEL_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/STA/FMC_C_W1.VEL'
PLT_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FMC_C_W1.PLT'

ALPSAP_SNX = readSNX(SNX_file, 'all')
% ALPSAP_CRD = readCRD(CRD_file);
% ALPSAP_VEL = readVEL(VEL_file);

%%
% flags = ALPSAP_CRD(:,7);
% range_flags = 1:length(flags);
% range_flags = range_flags(ismember(flags, {'A','W'}) == 1);  % not needed anymore

CRD_all =  ALPSAP_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL_all =  ALPSAP_SNX.SOLUTION.ESTIMATE.Data.VEL;
names_all = ALPSAP_SNX.SITE.ID.CODE;
names_all = cellstr(names_all);

% transform
[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);

SNX_cov = ALPSAP_SNX.SOLUTION.COVA_ESTIMATE;
% transform covariances from xyz to enu
[CovVenuSNX, SigmaVenu, CorrVen, AngleV] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'VEL');
[CovRenuSNX, SigmaRenu, CorrRen, AngleR] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'CRD');

%% remove Eurasia plate motion
% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Omega_Eur = [55.7892, -97.8099,   2.629e-07 ]';

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);
[V_res_xyz_all] = remove_plate_motion(CRD_all, VEL_all, Omega_Eur);
[Ve_res_all, Vn_res_all, Vu_res_all] = XYZ2ENU(CRD_all,V_res_xyz_all); % NEU components, [m/yr m/yr m/yr]

%% Megre artificial stations

[CRD,VEL,names] = merge_stations(CRD_all,VEL_all,names_all);
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);

[Ve,     Vn,     Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
[Ve_res, Vn_res, Vu_res]             = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

%% scale factor for presicion

SiteDome = [ALPSAP_SNX.SITE.ID.CODE, repmat(char(' '),size(CRD_all,1),1), ALPSAP_SNX.SITE.ID.DOMES];
SiteDome_list = cellstr(SiteDome);
[rmsENU(:,2), rmsENU(:,1), rmsENU(:,3)] = get_PLT_residuals(PLT_file, SiteDome_list);
%%
scalePLT_SNX_sig = mean(rmsENU,1) ./ mean(SigmaRenu,1);
disp('mean rms ,                  E         N         U')
disp(['PLT/SNX_sig           : ', num2str(scalePLT_SNX_sig,'%10.2f')])
mean(scalePLT_SNX_sig )

%% save  full horizontal velocity field
clc
% scaleFactor = 50*[1 1 1];
scaleFactor = scalePLT_SNX_sig/1;
VelocityField_hor = [long_all, lat_all, Ve_res_all, Vn_res_all, SigmaVenu(:,1)*scaleFactor(1), SigmaVenu(:,2)*scaleFactor(2),CorrVen];
writeVelocityFieldwithCovGMT(VelocityField_hor, names_all, '~/Alpen_Check/MAP/Velocity_field_horizontal.txt')
writeVelocityFieldVertical(long_all, lat_all, Vu_res_all, names_all, '~/Alpen_Check/MAP/Velocity_field_vertical.txt')

%% save Error Bars for GMT
fileID = fopen('~/Alpen_Check/MAP/CombNet/Vu_error_bars_all.txt', 'w');
fprintf(fileID, '# Velocity Field Error Bars Lat=Lat+Vu*scale, SigmaVu (mm/yr) -> SigmaVu[deg/yr] (for ploting with "gmt psxy -Ex" ) \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],      Sigma U [deg/yr],      \n');
formatStr = '%12.7f  %12.7f   %15e \n';
% d = diag(CovVenu);
SigmaVu = SigmaVenu(:,3) * scaleFactor(3); % scaled to abequate value [m/yr]
SigmaVu_deg_yr = SigmaVu*1000 * 0.36; 

%data = [long, lat + Vu_res * 1000 * 0.24, SigmaVu_deg_yr]; old

data = [long_all, lat_all + Vu_res_all * 1000 * 0.29, SigmaVu_deg_yr];

for i = 1:length(Vu_res_all)
   fprintf(fileID, formatStr, data(i,:)); 
end
fclose(fileID);
disp('Done')

%% run LSC for Combined network (ALP_NET + SAPOS)
%% for Horizontal 
clc
Outliers   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', ...
              'ROHR','BIWI','BI2I'};
iiOut = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);

V_enu_res = [Ve_res(iiSel), Vn_res(iiSel), Vu_res(iiSel)];
CovVenu2 = extractCovariance(CovVenuSNX, iiSel, [1 2 3], 'no split');
Cov_scale = 1;

clear LongGrid LatGrid V_def_tr1 rmsFit V_Sig_tr

%                                                                                                                           LongLim LatLing step dMax nMin flags ... 
[LongGrid, LatGrid, V_def_tr1, rmsFit, V_Sig_tr] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, Cov_scale, [0 18], [41 52], 0.5, 150, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'no filter');

%%
try
    close(fig1)
end
fig1 = figure(1);
% subplot(1,2,2)
hold on
sc=300
Earth_coast(2)
% text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all, Ve_res_all*sc,    Vn_res_all*sc,     0, 'b');
quiver(LongGrid,            LatGrid, V_def_tr1(:,1)*sc,V_def_tr1(:,2)*sc, 0, 'r', 'lineWidth',0.5)
xlim([0 19])
ylim([41 54])
title('Horizontal')


%%  Vertical
clc
Outliers = {'ELMO','WIEN','FERR', ...
            'BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE2','WT21','HKBL','PATK','PAT2', ...
            'HRIE','KTZ2', 'WLBH','ENTZ'};
        
iiOut    = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);
CovVenu2 = extractCovariance(CovVenuSNX, iiSel, [1 2 3], 'no split');
% Cov_scale = 5;
% VerrMin = (0.3/1000/Cov_scale)^2/1.8; % in [m/yr]^2, scaled to SNX
% CovVenu2(CovVenu2 < VerrMin) = VerrMin;
V_enu_res = [Ve_res_all(iiSel), Vn_res_all(iiSel), Vu_res_all(iiSel)];

% run for trend on regular grid

%                                                                                                                  LongLim LatLing step dMax nMin flags ... 
[LongGridv, LatGridv, V_def_v] = run_Collocation(long_all(iiSel), lat_all(iiSel), V_enu_res, CovVenu2, Cov_scale, [-1 18], [41 55], 0.5, 150, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%%
try
    close(fig2)
end
fig2 = figure(2);
% subplot(1,2,1)
hold on

Earth_coast(2)
text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all,  zeros(size(Vu_res_all)),    Vu_res_all*sc,   0, 'b');
quiver(LongGridv,           LatGridv, zeros(size(V_def_v(:,1))),  V_def_v(:,3)*sc, 0, 'r', 'lineWidth',0.5)
xlim([0 19])
ylim([41 52])
title('Vertical')

%% save results
writeDeformationFieldGMT([LongGrid,  LatGrid,  V_def_tr1(:,1),           V_def_tr1(:,2)], '~/Alpen_Check/MAP/CombNet/Def_field_Vh.txt','noCov')
writeDeformationFieldGMT([LongGridv, LatGridv, zeros(size(V_def_v(:,1))),V_def_v(:,3)],   '~/Alpen_Check/MAP/CombNet/Def_field_Vu.txt','noCov')

%% Run Kriging (to smooth LSC results for vertical component)
c = [-3:0.5:3];
[LongK_Stack, LatK_Stack, VuK_Stack, fig1, figx2, V_p] = runKrigingAtPoints(LongGridv, LatGridv, V_def_v, long ,lat, Vu_res, iiSel,c,5,[],0.1, names);

%% write results of kriging
write_xyzTable([LongK_Stack, LatK_Stack, VuK_Stack],    '~/Alpen_Check/MAP/CombNet/Vel_up_Kriging.txt', '%8.3f %8.3f %8e\n');

%% compute Strain 
DeformationField = [LongGrid,  LatGrid,  V_def_tr1(:,1:2)];
Strain = getStrainMap(DeformationField);

plotStrainNormal(Strain, 10^7*1);

%% write files with strain results

writeStrain2GMT(     Strain, '~/Alpen_Check/MAP/CombNet/Strain/StrainField.txt')
writeStrainSum2GMT(  Strain, '~/Alpen_Check/MAP/CombNet/Strain/StrainSum.txt')
writeStrainShear2GMT(Strain, '~/Alpen_Check/MAP/CombNet/Strain/StrainShear.txt')



