% combined multiyear ALP_NET + SAPOS = ALP_SAP
% from weekly solutions

%% SAPOS network
%
% 
SNX_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/SOL/FMC_C_W1.SNX'
CRD_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/STA/FMC_C_W1.CRD'
VEL_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/STA/FMC_C_W1.VEL'
PLT_file='/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FMC_C_W1.PLT'

ALPSAP_SNX = readSNX(SNX_file, 'all')
ALPSAP_CRD = readCRD(CRD_file);
ALPSAP_VEL = readVEL(VEL_file);

%%
flags = ALPSAP_CRD(:,7);
range_flags = 1:length(flags);
range_flags = range_flags(ismember(flags, {'A','W'}) == 1);


CRD_all =  ALPSAP_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL_all =  ALPSAP_SNX.SOLUTION.ESTIMATE.Data.VEL;
names_all = ALPSAP_SNX.SITE.ID.CODE;

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);

SNX_cov = ALPSAP_SNX.SOLUTION.COVA_ESTIMATE;
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

%% scale factor for presicion

SiteDome = [ALPSAP_SNX.SITE.ID.CODE, repmat(char(' '),length(range_flags),1), ALPSAP_SNX.SITE.ID.DOMES];
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
VelocityField = [long_all, lat_all, Ve_res_all, Vn_res_all, SigmaVenu(:,1)*scaleFactor(1), SigmaVenu(:,2)*scaleFactor(2),CorrVen];
writeVelocityFieldwithCovGMT(VelocityField, names_all, '/home/gast/Alpen_Check/MAP/Velocity_field_horizontal.txt')

%% save Error Bars for GMT
fileID = fopen('~/Alpen_Check/MAP/VelocityField/Vu_bars_all2.txt', 'w');
fprintf(fileID, '# Velocity Field Error Bars Lat=Lat+Vu*scale, SigmaVu (mm/yr) -> SigmaVu[deg/yr] (for ploting with "gmt psxy -Ex" ) \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],      Sigma U [deg/yr],      \n');
formatStr = '%12.7f  %12.7f   %15e \n';
% d = diag(CovVenu);
SigmaVu = SigmaVenu(:,3) * scaleFactor(3); % scaled to abequate value [m/yr]
SigmaVu_deg_yr = SigmaVu*1000 * 0.7; 

%data = [long, lat + Vu_res * 1000 * 0.24, SigmaVu_deg_yr]; old

data = [long_all, lat_all + Vu_res_all * 1000 * 0.46, SigmaVu_deg_yr];

for i = 1:length(Vu_res_all)
   fprintf(fileID, formatStr, data(i,:)); 
end
fclose(fileID);
disp('Done')

%% SAPOS Combined

clc
iiSel=[1:length(long_all)];
V_enu_res = [Ve_res_all(iiSel), Vn_res_all(iiSel), Vu_res_all(iiSel)];

CovVenu2 = extractCovariance(CovVenuSNX, iiSel, [1 2 3], 'no split');


% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD_all,VEL_all,names_all);


%% for Horizontal 
clc
Outliers   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', ...
              'ROHR','BIWI','BI2I'};
iiOut = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);
V_enu_res = [Ve_res_all(iiSel), Vn_res_all(iiSel), Vu_res_all(iiSel)];
CovVenu2 = extractCovariance(CovVenuSNX, iiSel, [1 2 3], 'no split');
Cov_scale = 20;
clear LongGrid LatGrid V_def_tr1 rmsFit V_Sig_tr
[LongGrid, LatGrid, V_def_tr1, rmsFit, V_Sig_tr] = run_Collocation(long_all(iiSel), lat_all(iiSel), V_enu_res, CovVenu2, Cov_scale, [1 16], [42 51], 0.25, 150, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%%
try
    close(fig1)
end
fig1 = figure(1);
% subplot(1,2,2)
hold on
Earth_coast(2)
% text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all, Ve_res_all*sc,      Vn_res_all*sc,    0, 'b');
quiver(LongGrid,            LatGrid, V_def_tr1(:,1)*sc,V_def_tr1(:,2)*sc,0, 'r', 'lineWidth',0.5)
xlim([0 17])
ylim([41 52])
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

%% run for trend on regular grid
[LongGridv, LatGridv, V_def_v] = run_Collocation(long_all(iiSel), lat_all(iiSel), V_enu_res, CovVenu2, Cov_scale, [0 18], [42 53], 0.5, 200, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%%
try
    close(fig2)
end
fig2 = figure(2);
% subplot(1,2,1)
hold on
sc=300
Earth_coast(2)
% text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all,  zeros(size(Vu_res_all)),    Vu_res_all*sc,   0, 'b');
quiver(LongGridv,           LatGridv, zeros(size(V_def_v(:,1))),  V_def_v(:,3)*sc, 0, 'r', 'lineWidth',0.5)
xlim([0 17])
ylim([41 52])
title('Vertical')







