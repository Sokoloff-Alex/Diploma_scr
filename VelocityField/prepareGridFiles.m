% prepare grid field 

% FMC_IGB_W8_SNX = readSNX('../dat/SNX/FMC_IGB_W8.SNX','noCov');
% save('../dat/SNX/FMC_IGB_W8_SNX.mat','FMC_IGB_W8_SNX');

load('../dat/SNX/FMC_IGB_W8_SNX.mat')

%% write grid 

x = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD(:,1);
y = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD(:,2);
z = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD(:,3);

s_x = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD_STD(:,1);
s_y = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD_STD(:,2);
s_z = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.CRD_STD(:,3);

vx = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL(:,1);
vy = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL(:,2);
vz = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL(:,3);

s_vx = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL_STD(:,1);
s_vy = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL_STD(:,2);
s_vz = FMC_IGB_W8_SNX.SOLUTION.ESTIMATE.Data.VEL_STD(:,3);

CODE  = FMC_IGB_W8_SNX.SOLUTION.EPOCHS.CODE;
PT    = FMC_IGB_W8_SNX.SOLUTION.EPOCHS.PT;
SolN  = FMC_IGB_W8_SNX.SOLUTION.EPOCHS.SolN;
Start = FMC_IGB_W8_SNX.SOLUTION.EPOCHS.DATA_START;
End   = FMC_IGB_W8_SNX.SOLUTION.EPOCHS.DATA_END;

DomesN= FMC_IGB_W8_SNX.SITE.ID.DOMES;

%% approx duration

t_start = [str2num(Start(:,1:2)) + str2num(Start(:,4:6))/365.25]
t_end   = [str2num(End(  :,1:2)) + str2num(End(  :,4:6))/365.25]

dt = t_end - t_start;

%% CRD_XYZ
clc
SigmaScale = [30 30 30];
data = [x s_x*SigmaScale(1) y s_y*SigmaScale(2) z s_z*SigmaScale(3)];
formatStr  = '%3d  %4s %9s  %15.5f %8.5f  %15.5f %8.5f  %15.5f %8.5f   %2s   %2d    20%6s  20%6s\n';
disp('Nr.  STATION NAME            X[m]     sig_X[m]         Y[m]     sig_Y[m]         Z[m]     sig_Z[m]   ID-SNX     START      END')
for i = 1:length(x)
    fprintf(formatStr, i, CODE(i,:), DomesN(i,:), data(i,:), PT(i,:), SolN(i,:), Start(i,1:6), End(i,1:6)); 
end

%% VEL_XYZ
clc
SigmaScale = [30 30 30];
data = [vx s_vx*SigmaScale(1) vy s_vy*SigmaScale(2) vz s_vz*SigmaScale(3)];
formatStr  = '%3d  %4s %9s  %10.5f %8.5f   %10.5f %8.5f   %10.5f %8.5f     %2s   %2d    20%6s  20%6s\n';
disp('Nr.  STATION NAME       VX[m/a]  sig_VX[m/a]  VY[m/a]  sig_VY[m/a]  VZ[/am]  sig_VZ[m/a]  ID-SNX     START      END')
for i = 1:length(x)
    fprintf(formatStr, i, CODE(i,:), DomesN(i,:), data(i,:), PT(i,:), SolN(i,:), Start(i,1:6), End(i,1:6)); 
end

%%
[Ve,Vn, Vu, lat,lon, h] = XYZ2ENU([x y z], [vx vy vz]);
[sVe,sVn, sVu] = XYZ2ENU([x y z], [s_vx s_vy s_vz]);
[se,sn, su] = XYZ2ENU([x y z], [s_x s_y s_z]);

%% CRD_XYZ
clc
SigmaScale = [44 28 19];
data = [lon se*SigmaScale(1) lat sn*SigmaScale(2) h su*SigmaScale(3)];
formatStr  = '%3d  %4s %9s  %10.7f %8.5f     %10.7f  %8.5f   %10.5f %8.5f     %2s   %2d    20%6s  20%6s\n';
disp('Nr.  STATION NAME    Long[deg]  sig_Long[m]  Lat[deg]    sig_Lat[m]  Height[m]  sig_H[m]     ID-SNX     START      END')
for i = 1:5 %length(x)
    fprintf(formatStr, i, CODE(i,:), DomesN(i,:), data(i,:), PT(i,:), SolN(i,:), Start(i,1:6), End(i,1:6)); 
end

%% VEL_XYZ
clc
SigmaScale = [44 28 19];
data = [Ve sVe*SigmaScale(1) Vn sVn*SigmaScale(2) Vu sVu*SigmaScale(3)];
formatStr  = '%3d  %4s %9s  %10.5f %8.5f   %10.5f %8.5f   %10.5f %8.5f     %2s   %2d    20%6s  20%6s\n';
disp('Nr.  STATION NAME       Ve[m/a]  sig_Ve[m/a]  Vn[m/a]  sig_Vn[m/a]  Vu[/am]  sig_Vu[m/a]  ID-SNX     START      END')
for i = 1:5 %length(x)
    fprintf(formatStr, i, CODE(i,:), DomesN(i,:), data(i,:), PT(i,:), SolN(i,:), Start(i,1:6), End(i,1:6)); 
end



