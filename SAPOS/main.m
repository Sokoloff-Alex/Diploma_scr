%% SAPOS network
%
% 
SNX_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/SOL/FMC_SAP_1.SNX'
CRD_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FMC_SAP_1.CRD'
VEL_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FMC_SAP_1.VEL'
PLT_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/OUT/FMC_SAP_1.PLT'

SAPOS_SNX = readSNX(SNX_file, 'all')
SAPOS_CRD = readCRD(CRD_file);
SAPOS_VEL = readVEL(VEL_file);


%
flags = SAPOS_CRD(:,7);
range_flags = 1:length(flags);
range_flags = range_flags(ismember(flags, {'A','W'}) == 1);


%CRD_all = cell2mat( SAPOS_CRD(range_flags,4:6));
%VEL_all = cell2mat( SAPOS_VEL(range_flags,4:6));
%names_all = SAPOS_CRD(range_flags,2);
%DOMES = SAPOS_CRD(range_flag,3);


CRD_all =  SAPOS_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL_all =  SAPOS_SNX.SOLUTION.ESTIMATE.Data.VEL;
names_all = SAPOS_SNX.SITE.ID.CODE;


[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);


SNX_cov = SAPOS_SNX.SOLUTION.COVA_ESTIMATE;
[CovVenuSNX, SigmaVenu, CorrVen, AngleV] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'VEL');
[CovRenuSNX, SigmaRenu, CorrRen, AngleR] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'CRD');

% remove Eurasia plate motion
% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Omega_Eur = [55.7892, -97.8099,   2.629e-07 ]';

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);
[V_res_xyz_all] = remove_plate_motion(CRD_all, VEL_all, Omega_Eur);
[Ve_res_all, Vn_res_all, Vu_res_all] = XYZ2ENU(CRD_all,V_res_xyz_all); % NEU components, [m/yr m/yr m/yr]

%% scale factor for presicion

SiteDome = [SAPOS_SNX.SITE.ID.CODE, repmat(char(' '),104,1), SAPOS_SNX.SITE.ID.DOMES];
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
scaleFactor = scalePLT_SNX_sig/1.5;
VelocityField = [long_all, lat_all, Ve_res_all, Vn_res_all, SigmaVenu(:,1)*scaleFactor(1), SigmaVenu(:,2)*scaleFactor(2),CorrVen];
writeVelocityFieldwithCovGMT(VelocityField, names_all, '/home/gast/SAPOS/MAP/Velocity_field_horizontal.txt')

%% save Error Bars for GMT
fileID = fopen('~/SAPOS/MAP/VelocityField/Vu_bars_all2.txt', 'w');
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

