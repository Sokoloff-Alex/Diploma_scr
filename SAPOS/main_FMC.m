% main_FNC

SNX_file_fnc='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/SOL/FNC_SAP_1.SNX'
CRD_file_fnc='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FNC_SAP_1.CRD'
VEL_file_fnc='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/STA/FNC_SAP_1.VEL'
PLT_file_fnc='/home/gast/GPSDATA/CAMPAIGN52/SAPOS/OUT/FNC_SAP_1.PLT'

SAPOS_fnc_SNX = readSNX(SNX_file_fnc, 'all')
SAPOS_fnc_CRD = readCRD(CRD_file_fnc);
SAPOS_fnc_VEL = readVEL(VEL_file_fnc);

%%
flags_fnc = SAPOS_fnc_CRD(:,7);
range_flags_fnc = 1:length(flags_fnc);
range_flags_fnc = range_flags_fnc(ismember(flags_fnc, {'A','W'}) == 1);


%CRD_all = cell2mat( SAPOS_CRD(range_flags,4:6));
%VEL_all = cell2mat( SAPOS_VEL(range_flags,4:6));
%names_all = SAPOS_CRD(range_flags,2);
%DOMES = SAPOS_CRD(range_flag,3);


CRD_all_fnc   = SAPOS_fnc_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL_all_fnc   = SAPOS_fnc_SNX.SOLUTION.ESTIMATE.Data.VEL;
names_all_fnc = SAPOS_fnc_SNX.SITE.ID.CODE;


[Ve_all_fnc,Vn_all_fnc, Vu_all_fnc, lat_all_fnc, long_all_fnc,  h_all_fnc]  = XYZ2ENU(CRD_all_fnc,VEL_all_fnc);


% remove Eurasia plate motion
% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Omega_Eur = [55.7892, -97.8099,   2.629e-07 ]';

[Ve_all_fnc,Vn_all_fnc, Vu_all_fnc, lat_all_fnc, long_all_fnc,  h_all_fnc]  = XYZ2ENU(CRD_all_fnc,VEL_all_fnc);
[V_res_xyz_all_fnc] = remove_plate_motion(CRD_all_fnc, VEL_all_fnc, Omega_Eur);
[Ve_res_all_fnc, Vn_res_all_fnc, Vu_res_all_fnc] = XYZ2ENU(CRD_all_fnc,V_res_xyz_all_fnc); % NEU components, [m/yr m/yr m/yr]


%% exclude '1277'
Ve_all_fnc = Ve_all_fnc([1:62,64:end]);
Vn_all_fnc = Vn_all_fnc([1:62,64:end]);
Vu_all_fnc = Vu_all_fnc([1:62,64:end]);
Ve_res_all_fnc = Ve_res_all_fnc([1:62,64:end]);
Vn_res_all_fnc = Vn_res_all_fnc([1:62,64:end]);
Vu_res_all_fnc = Vu_res_all_fnc([1:62,64:end]);
lat_all_fnc = lat_all_fnc([1:62,64:end]);
long_all_fnc = long_all_fnc([1:62,64:end]);
names_all_fnc = names_all_fnc([1:62,64:end],:);

%%
try
    close(fig1)
end
fig1 = figure(1);
hold on
sc=1000
% text(wrapTo180(long_all_fnc), lat_all_fnc, names_all_fnc)
bias_Vup = mean(Vu_res_all_fnc)- 0.0006;
dvup = Vu_res_all_fnc - bias_Vup - Vu_res_all;
% quiver(wrapTo180(long_all_fnc), lat_all_fnc, zeros(size(Vu_res_all_fnc)), (Vu_res_all_fnc - bias_Vup)*sc, 0, 'b')
% quiver(wrapTo180(long_all_fnc), lat_all_fnc, zeros(size(Vu_res_all_fnc)), (Vu_res_all)*sc,0,'Color',[200 0 200]/256)
quiver(wrapTo180(long_all_fnc), lat_all_fnc, zeros(size(Vu_res_all_fnc)), dvup*sc,0,'r')

%%
% close all
try
    close(fig2);
end
fig2 = figure(2);
hold on
sc=1000;
% text(wrapTo180(long_all_fnc), lat_all_fnc, names_all_fnc)
bias_Ve = mean(Ve_res_all_fnc);
bias_Vn = mean(Vn_res_all_fnc);
dve = Ve_res_all_fnc - bias_Ve - Ve_res_all;
dvn = Vn_res_all_fnc - bias_Vn - Vn_res_all;

% quiver(wrapTo180(long_all_fnc), lat_all_fnc, (Ve_res_all_fnc - bias_Ve)*sc, (Vn_res_all_fnc - bias_Vn)*sc, 0, 'b')
% quiver(wrapTo180(long_all_fnc), lat_all_fnc, (Ve_res_all)*sc, (Vn_res_all)*sc,0,'Color',[200 0 200]/256)
quiver(wrapTo180(long_all_fnc), lat_all_fnc, dve*sc, dvn*sc,0,'r')

