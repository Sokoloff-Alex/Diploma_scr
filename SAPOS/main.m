%% SAPOS network
%
% 
SNX_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/SOL/FMC_SAP_1.SNX'
CRD_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FMC_SAP_1.CRD'
VEL_file='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FMC_SAP_1.VEL'

SAPOS_SNX = readSNX(SNX_file, 'all')
SAPOS_CRD = readCRD(CRD_file);
SAPOS_VEL = readVEL(VEL_file);


%%
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

%% remove Eurasia plate motion
% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Omega_Eur = [55.7892, -97.8099,   2.629e-07 ]';

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);
[V_res_xyz_all] = remove_plate_motion(CRD_all, VEL_all, Omega_Eur);
[Ve_res_all, Vn_res_all, Vu_res_all] = XYZ2ENU(CRD_all,V_res_xyz_all); % NEU components, [m/yr m/yr m/yr]

%% save  full horizontal velocity field
clc
scaleFactor = 20*[1 1 1];
VelocityField = [long_all, lat_all, Ve_res_all, Vn_res_all, SigmaVenu(:,1)*scaleFactor(1), SigmaVenu(:,2)*scaleFactor(2),CorrVen];
writeVelocityFieldwithCovGMT(VelocityField, names_all, '/home/gast/SAPOS/MAP/Velocity_field_horizontal.txt')

