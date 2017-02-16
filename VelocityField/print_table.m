% print table

close all
clear all
clc

%% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Orogen_Alp  = importOrogen('dat/PB2002_orogen_Alps.txt');
lwmask = struct2array(load('dat/lwmask25.mat'));

% [Etopo_Europe, refvec_Etopo] = etopo('../../../../MAP/etopo1_bed_c_f4/etopo1_bed_c_f4.flt', 1, [40 54], [-7 19]); % ETOPO
load('ETOPO_Alps.mat')
% Adriatics = struct2array(load('dat/Adriatics.mat'));
%
ALP_NET_CRD = readCRD('STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('STA/FMC_IGB_W7.VEL');


%
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

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);

% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD_all,VEL_all,names_all);

% remove Eurasia plate motion
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);
[V_res_xyz_all] = remove_plate_motion(CRD_all, VEL_all, Omega_Eur);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]
[Ve_res_a, Vn_res_a, Vu_res_a] = XYZ2ENU(CRD_all,V_res_xyz_all); % NEU components, [m/yr m/yr m/yr]


%% Covariance

load('../dat/SNX/SINEX.mat')
SNX_cov = SINEX.SOLUTION.COVA_ESTIMATE;
[CovVenuSNX, SigmaVenu, CorrVen, AngleV] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'VEL');
[CovRenuSNX, SigmaRenu, CorrRen, AngleR] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'CRD');
[CovVenu] = megreCov(CovVenuSNX, names_all);
[CRD, SigmaVenu_merged, name ] = merge_stations(CRD_all,SigmaVenu,names_all);
[CRD, AngleV_merged] = merge_stations(CRD_all,AngleV,names_all);
Angle_v = AngleV_merged(:,1);

DOMES_all = SINEX.SITE.ID.DOMES;


%% compute common observation period
% 
t_start = SINEX.SOLUTION.EPOCHS.DATA_START;
t_end   = SINEX.SOLUTION.EPOCHS.DATA_END;

t_start = [str2num(t_start(:,1:2)) + str2num(t_start(:,4:6))/365.25]
t_end   = [str2num(t_end(  :,1:2)) + str2num(t_end(  :,4:6))/365.25]

dt = t_end - t_start;
[dtobs_sum] = merge_stations_sum(dt,names_all);

%%
[t_obs_start] = merge_stations_param(t_start, names_all, 'min');
[t_obs_end  ] = merge_stations_param(t_end,   names_all, 'max');
[CorrVen_ave, names_1, names_2] = merge_stations_param(CorrVen, names_all, 'mean');

Sigma_Venu = SigmaVenu_merged * sqrt(1.8) * 20 * 1000; % im mm/yr


%% save full table of TOTAL velocity field in compact format

clc
dt2 = (2000+t_obs_end) - (2000+t_obs_start);
td_err = dt2 - dtobs_sum;

data = [wrapTo180(long), lat, h, Ve*1000, Sigma_Venu(:,1), Vn*1000, Sigma_Venu(:,2), Ve_res*1000, Sigma_Venu(:,1), Vn_res*1000, Sigma_Venu(:,2), Vu*1000, Sigma_Venu(:,3), 2000+t_obs_start, 2000+t_obs_end,  dtobs_sum];
formatStr = '%4s  %8.2f    %8.2f     %6.1f        %5.2f ± %4.2f     %5.2f ± %4.2f     %5.2f ± %4.2f     %5.2f ± %4.2f     %5.2f ± %4.2f    %5.1f - %5.1f  %4.1f \n';
disp('Site     Long_app,   Lat_app,  Height_app     Ve      SVe      Vn      SVn      Ve_res  SVe      Vn_res  SVn      Vu     Svu      Start    End    Tobs  ')
for i = 1:length(long)
    fprintf(formatStr, name{i}, data(i,:)); 
end

%% save full table of RESIDUAL velocity field in compact format

clc
dt2 = (2000+t_obs_end) - (2000+t_obs_start);
td_err = dt2 - dtobs_sum;
data = [wrapTo180(long), lat, h, Ve_res*1000, Sigma_Venu(:,1), Vn_res*1000, Sigma_Venu(:,2), Vu_res*1000, Sigma_Venu(:,3), 2000+t_obs_start, 2000+t_obs_end,  dtobs_sum];
formatStr = '%4s  %8.2f    %8.2f     %6.1f        %5.2f ± %4.2f     %5.2f ± %4.2f     %5.2f ± %4.2f    %5.1f - %5.1f  %4.1f \n';
disp('Site     Long_app,   Lat_app,  Height_app     Ve_res  SVe      Vn_res  SVn      Vu     Svu      Start    End    Tobs  ')
for i = 1:length(long)
    fprintf(formatStr, name{i}, data(i,:)); 
end

%% save full table of velocity field in extended format

clc
data = [wrapTo180(long), lat, h, Ve_res*1000, Sigma_Venu(:,1), Vn_res*1000, Sigma_Venu(:,2), Vu_res*1000, Sigma_Venu(:,3), 2000+t_obs_start, 2000+t_obs_end,  dtobs_sum];
formatStr = '%4s  %14.7f    %14.7f     %12.5f        %5.2f ± %4.2f     %5.2f ± %4.2f     %5.2f ± %4.2f    %5.1f - %5.1f  %4.1f \n';
disp('Site     Long_app,   Lat_app,  Height_app     Ve_res  SVe      Vn_res  SVn      Vu     Svu      Start    End    Tobs  ')
for i = 1:length(long)
    fprintf(formatStr, name{i}, data(i,:)); 
end

%%
clc
n  = length(lat_all);
SigmaVenu_all = SigmaVenu * sqrt(1.8) * 20 * 1000; 
data = [wrapTo180(long_all), lat_all, h_all, Ve_res_a*1000, SigmaVenu_all(:,1), Vn_res_a*1000, SigmaVenu_all(:,2), Vu_res_a*1000, SigmaVenu_all(:,3), 2000+t_start, 2000+t_end,  dt];
formatStr = '%4s  %4s %9s %7.2f %8.2f  %6.1f  %5.2f ± %4.2f  %5.2f ± %4.2f  %5.2f ± %4.2f    %5.1f - %5.1f  %4.1f \n';
disp('SITE: Site__________       Longitude        Latitude           Height           Ve_res  SVe      Vn_res  SVn      Vu     Svu      Start    End    Tobs  ')
for i = 1:length(long_all)
    fprintf(formatStr, names_2{i},  names_all{i},DOMES_all(i,:) , data(i,:)); 
    if i<n
        if (~strcmp(names_2{i}, names_2{i+1})) 
            disp('--------------------')
        end
    end

end





