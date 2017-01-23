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

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);

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
[CRD, AngleV_merged] = merge_stations(CRD_all,AngleV,names_all);
Angle_v = AngleV_merged(:,1);

%% 

%% save Error Bars for GMT
fileID = fopen('~/Alpen_Check/MAP/VelocityField/Vu_bars_all.txt', 'w');
fprintf(fileID, '# Velocity Field Error Bars Lat=Lat+Vu*scale, SigmaVu (mm/yr) -> SigmaVu[deg/yr] (for ploting with "gmt psxy -Ex" ) \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],      Sigma U [deg/yr],      \n');
formatStr = '%12.7f  %12.7f   %15e \n';
d = diag(CovVenu);
SigmaVu = SigmaVenu_merged(:,3) * 1.8^(1/2) * 20; % scaled to abequate value [m/yr]
SigmaVu_deg_yr = SigmaVu*1000 * 0.36; 

%data = [long, lat + Vu_res * 1000 * 0.24, SigmaVu_deg_yr]; old

data = [long, lat + Vu_res * 1000 * 0.29, SigmaVu_deg_yr];

for i = 1:length(Vu_res)
   fprintf(fileID, formatStr, data(i,:)); 
end
fclose(fileID);

%% compute common observation period
% 
% t_start = SINEX.SOLUTION.EPOCHS.DATA_START;
% t_end   = SINEX.SOLUTION.EPOCHS.DATA_END;
% 
% t_start = [str2num(t_start(:,1:2)) + str2num(t_start(:,4:6))/365.25]
% t_end   = [str2num(t_end(  :,1:2)) + str2num(t_end(  :,4:6))/365.25]
% 
% dt = t_end - t_start;
% 
% %%
% [dtobs_sum] = merge_stations_sum(dt,names_all);

%%
% close all
% figure(1)
% hist(dtobs_sum,10)

%% save full table of velocity field
% % d = diag(CovVenu);
% clc
% Sigma_Venu = SigmaVenu_merged * 1.8^(1/2) * 20 * 1000;
% data = [wrapTo180(long), lat, Ve_res*1000, Sigma_Venu(:,1), Vn_res*1000, Sigma_Venu(:,2), Angle_v, Vu_res*1000, Sigma_Venu(:,3), dtobs_sum];
% formatStr = '%4s  %8.3f\t %8.3f\t %5.2f  %5.2f   %5.2f  %5.2f   %5.1f     %5.2f  %5.2f   %4.1f \n';
% disp('Site  Long,     Lat,          Ve     SVe     Vn     SVn     angle     Vu     Svu   Tobs')
% for i = 1:length(long)
%     fprintf(formatStr, name{i}, data(i,:)); 
% end


%% save  full horizontal velocity field
clc
fileID = fopen('~/Alpen_Check/MAP/VelocityField/Vel_hor_total_SNX_nonames.txt','w');
Sigma_Venu = SigmaVenu_merged * 1.8^(1/2) * 20 * 1000;
data = [wrapTo180(long_all), lat_all, Ve_all, Vn_all, SigmaVenu(:,1)*1.8*20, SigmaVenu(:,2)*1.8*20, AngleV];
formatStr = '%10.7f  %10.7f   %8.5f  %8.5f   %10.6f  %10.6f  %5.1f \n';
header = 'Long [deg],  Lat[deg],    Ve[m/yr],  Vn[m/yr],  SVn[m/yr], SVe[m/yr], angle[deg],    Name\n';
fprintf(fileID, header);
for i = 1:length(names_all)
    fprintf(fileID, formatStr, data(i,:)); 
end
fclose(fileID);
%% for Horizontal 

Outliers   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', ...
                'ROHR','BIWI','BI2I'};
iiOut = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);

%% wrie outliers files
% writeVelocityFieldVerticalOutliers2GMT(long(iiOut), lat(iiOut), Vu_res(iiOut), names(iiOut), '~/Alpen_Check/MAP/VelocityField/VelocityFieldVertical_out.txt');

%%
clc
CovVenu2 = extractCovariance(CovVenu, iiSel, [1 2 3], 'no split');
V_enu_res = [Ve_res(iiSel), Vn_res(iiSel), Vu_res(iiSel)];
%% for Horizontal
clc
%                                                                                                             LongLim LatLing step dMax nMin flags ...                                         
[LongGrid, LatGrid, V_def, rmsFit, V_SigPred] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [-4 18], [41 53], 1, 150, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'no filter');

%% for Vertical       

% Block selection for Vertical

Outliers   = {'HELM', 'WIEN', 'OGAG', 'OBE2', ...
              'ROHR','BIWI','BI2I', 'MANS'};
            
Outliers = {'ELMO','WIEN','FERR', ...
            'BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE2','WT21','HKBL','PATK','PAT2', ...
            'HRIE','KTZ2', 'WLBH'};
        
iiOut    = selectRange(names, Outliers);
iiSel = setdiff([1:198], iiOut);
CovVenu2 = extractCovariance(CovVenu, iiSel, [1 2 3], 'no split');
VerrMin = (0.3/1000/25)^2/1.8; % in [m/yr]^2, scaled to SNX
CovVenu2(CovVenu2 < VerrMin) = VerrMin;
V_enu_res = [Ve_res(iiSel), Vn_res(iiSel), Vu_res(iiSel)];
[LongGrid, LatGrid, V_def, rmsFit, V_SigPred, Cs0] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [-4 18], [41 53], 1, 250, 10, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%% run Kriging
close all
c = [-2:0.5:-.5, 0.5:0.5:3];
runKriging(LongGrid, LatGrid, V_def, long ,lat, Vu_res, iiSel,c,5);

%% run Kriging
close all
sig = sqrt(diag(CovVenu));
sig = [sig(1:3:end),sig(2:3:end),sig(1:3:end)]*sqrt(1.8)*50;
sig(sig < 0.0001) = 0.0001;
c = [0.2:0.2:2];
[LongK_Stack, LatK_Stack, VuK_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, V_SigPred/1000, long ,lat, sig(:,3), iiSel,c,5);
%%
write_xyzTable([LongK_Stack, LatK_Stack, VuK_Stack],'~/Alpen_Check/MAP/Deformation/Verr_up_est.txt','%12.7f %12.7f %12.3f\n')


%%
print(fig2,'-depsc','Verr_est.epsc','-r300')
print(fig2,'-dpng', 'Verr_est.png', '-r300')

%% check individual
run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [10 10], [45 45], 1, 250, 10, 'exp1', '-v', 'no bias', 'tail 0', 'no corr', 'filter','plot');

%% Plot Results 2D map
clc
s = 500;  % [mm/yr]
s2 = 0.2;
try
    close (fig7)
end
clr = lines(8);
fig7 = figure(7);
hold on
xlim([-2 18])
ylim([41 52])
etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
% Earth_coast(2)
% plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
% plot(Adriatics(:,1),Adriatics(:,2) , '--k')
quiver(long(iiSel), lat(iiSel), Ve_res(iiSel)*s,    Vn_res(iiSel)*s,  0, 'r', 'lineWidth',0.5)
% quiver(LongGrid,      LatGrid,      V_def(:,1)*s,   V_def(:,2)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% errorbar(long(iiSel), lat(iiSel)+ Vu_res(iiSel)*s, SigmaVenu(iiSel,3)*s, '.m')
% quiver(long(iiSel), lat(iiSel), zeros(size(iiSel))', Vu_res(iiSel)*s,  0, 'r', 'lineWidth',1)
quiver(long(iiOut),lat(iiOut), Ve_res(iiOut)*s,   Vn_res(iiOut)*s, 0, 'm', 'lineWidth',0.5)
% quiver(long(iiOut),lat(iiOut),zeros(size(iiOut))',  Vu_res(iiOut)*s, 0, 'm', 'lineWidth',1)
% text(long(iiOut),lat(iiOut), names(iiOut))
quiver(LongGrid,      LatGrid,      V_def(:,1)*s,   V_def(:,2)*s,    0, 'Color',clr(1,:), 'lineWidth',0.5)
% errorbar(LongGrid, LatGrid + V_def(:,3)*s, rmsFit(:,3)*s, '.k')
% errorbar(LongGrid, LatGrid + V_def(:,3)*s, V_SigPred(:,3)*s, '.k')
% quiver(LongGrid,      LatGrid,      zeros(size(LongGrid)),   V_def(:,3)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% text(long,lat, names)
plotErrorElipses('Cov', CovVenu*sqrt(s2),                         long,     lat,     Ve_res,     Vn_res,     s, 0.95, 'r')
plotErrorElipses('Sig', sqrt([V_SigPred(:,1),V_SigPred(:,2)])*s2, LongGrid, LatGrid, V_def(:,1), V_def(:,2), s, 0.95, 'b')
% plotErrorElipses('Sig', ([rmsFit(:,1),rmsFit(:,2)])*1000,  LongGrid, LatGrid, V_def(:,1), V_def(:,2), s, 0.95, 'm')


title('velocity field / deformation model') 
% xlabel('Velocity EW, [mm/yr]')
% ylabel('Velocity SN, [mm/yr]')
legend('Residual velocity','Outliers','LSC velocity','location','NorthWest')
set(gcf, 'PaperType', 'A4');
% print(fig7,'-dpng','-r300','Noise_Prop.png')


%% Plot for vertical
clc
close all
figure(2)
hold on

% Earth_coast(2)
% etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);

% Verr_grid = stack2grid([LongGrid, LatGrid, V_SigPred(:,3)]);
% quiver(long(iiSel), lat(iiSel), zeros(size(iiSel))',   Vu_res(iiSel)*s,  0, 'r')
% quiver(LongGrid,    LatGrid,    zeros(size(LongGrid)), V_def(:,3)*s,     0, 'Color',clr(1,:))

%
Cs0 = abs(Cs0);
V_SigPred = abs(V_SigPred);

a2 = abs(Cs0.^2 - V_SigPred.^2);
a = sqrt(abs(a2));
noise_assumed = 0.5 ; % [mm/yr]
Sig_conv = sqrt(noise_assumed^2 - a2);


subplot(2,2,1)
hold on
scatter(LongGrid, LatGrid, abs(Cs0(:,3))*200, (Cs0(:,3)),'fill')
legend('Cs0 only')
colorbar
% caxis([0 1.2])

subplot(2,2,2)
hold on
scatter(LongGrid, LatGrid, abs(V_SigPred(:,3))*200, (V_SigPred(:,3)),'fill')
legend('Predicted, LSC, Cs0 - a')
colorbar
% caxis([0 1])

subplot(2,2,3)
hold on
scatter(LongGrid, LatGrid, abs(a(:,3))*200, (a(:,3)),'fill','o')
legend('a only')
colorbar
% caxis([0 1.2])

subplot(2,2,4)
hold on
scatter(LongGrid, LatGrid, abs(rmsFit(:,3))*1000*200, rmsFit(:,3)*1000,'fill','o')
legend('sqrt(rms) fitting')
colorbar
% caxis([0 1])

%%
clc
close all
[LongK_Stack, LatK_Stack, Cs0_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, Cs0/1000, long ,lat, sig(:,3), iiSel,c,5);
% scatter(LatK_Stack, LongK_Stack, Cs0_Stack)
write_xyzTable([LongK_Stack, LatK_Stack, abs(Cs0_Stack)],'~/Alpen_Check/MAP/Deformation/Verr_up_Cs0.txt','%12.3f %12.3f %12e\n')

%%
close all
[LongK_Stack, LatK_Stack, V_SigPred_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, V_SigPred/1000, long ,lat, sig(:,3), iiSel,c,5);
write_xyzTable([LongK_Stack, LatK_Stack, abs(V_SigPred_Stack)],'~/Alpen_Check/MAP/Deformation/Verr_up_Mikhail.txt','%12.3f %12.3f %12e\n')

%%
close all
[LongK_Stack, LatK_Stack, Sig_conv_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, Sig_conv/1000, long ,lat, sig(:,3), iiSel,c,5);
% write_xyzTable([LongK_Stack, LatK_Stack, abs(Sig_conv_Stack)],'~/Alpen_Check/MAP/Deformation/Verr_up_conv.txt','%12.3f %12.3f %12e\n')

%%
close all
[LongK_Stack, LatK_Stack, rms_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, rmsFit, long ,lat, sig(:,3), iiSel,c,5);
write_xyzTable([LongK_Stack, LatK_Stack, abs(rms_Stack)],'~/Alpen_Check/MAP/Deformation/Verr_up_rms.txt','%12.3f %12.3f %12e\n')

%%
close all
[LongK_Stack, LatK_Stack, a_Stack , fig1, fig2] = runKriging(LongGrid, LatGrid, a/1000, long ,lat, sig(:,3), iiSel,c,5);
write_xyzTable([LongK_Stack, LatK_Stack, abs(a_Stack)],'~/Alpen_Check/MAP/Deformation/Verr_up_cov.txt','%12.3f %12.3f %12e\n')


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

