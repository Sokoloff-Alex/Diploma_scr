%%  Determine Scale factor 

close all
clear all
clc

%%
SINEX = readSNX('STA/FMC_IGB_W7.SNX','All');
% save('SINEX.mat','SINEX')
Lat = SINEX.SITE.ID.LAT;
Lon = SINEX.SITE.ID.LON;
CRD = SINEX.SOLUTION.ESTIMATE.Data.CRD;
CRD_std = SINEX.SOLUTION.ESTIMATE.Data.CRD_STD;
Sites_list = SINEX.SITE.ID.CODE;
Stations = cellstr(SINEX.SITE.ID.CODE);
DOMES    = cellstr(SINEX.SITE.ID.DOMES);
SiteDome = [SINEX.SITE.ID.CODE, repmat(char(' '),297,1), SINEX.SITE.ID.DOMES];
SiteDome_list = cellstr(SiteDome);

%%

clc
SiteDome_list_test = SiteDome_list(1:20);
[rmsENU(:,2), rmsENU(:,1), rmsENU(:,3), Sites] = get_PLT_residuals('../STA/FMC_IGB_W7.PLT', SiteDome_list);

% save('PLT_rmsENU.mat','rmsENU')

%%
find(strcmp(SiteDome_list,'KA3L 14216M001'))
range_w = [101 102 103 149 160 291 292 296];
Stations(range_w)

sel_good = [1 97 130 132 133 138 139 ]
Stations(sel_good)

%%
% wrong
% [CRD_std_enu(:,1), CRD_std_enu(:,2), CRD_std_enu(:,3)] = XYZ2ENU(CRD, CRD_std);
%% Transform sigma / covariance from XYZ to ENU
sigma_enu = zeros(size(CRD_std));
angle = size(CRD_std,1);
corr_en = size(CRD_std,1);
for i = 1:size(CRD_std,1)
    cov_xyz = [CRD_std(i,1)^2 0 0 ;
               0 CRD_std(i,2)^2 0 ;
               0 0 CRD_std(i,3)^2];
    [cov_enu, sigma_enu(i,:)] = covXYZ2ENU(cov_xyz, lat(1), long(1));
    corr_en(i,:) = cov_enu(1,2) / (sigma_enu(1) * sigma_enu(2));
    angle(i,:) = 90 -1/2 * atand(2*cov_enu(1,2) / (cov_enu(1) - cov_enu(2)));
end

CRD_std_enu = sigma_enu;

%%

SigmaRenu_Cov = SigmaRenu;
close all
try
    close (fig1)
end
fig1 = figure(1);
subplot(2,1,1)
hold on
plot(CRD_std*1000)
plot(SigmaRenu_Cov*1000*1.8,'.-')
plot([1, length(CRD_std_enu)], [mean(CRD_std_enu(:,1)), mean(CRD_std_enu(:,1))]*1000, '--b')
plot([1, length(CRD_std_enu)], [mean(CRD_std_enu(:,2)), mean(CRD_std_enu(:,2))]*1000, '--g')
plot([1, length(CRD_std_enu)], [mean(CRD_std_enu(:,3)), mean(CRD_std_enu(:,3))]*1000, '--r')
ylabel('[mm]')
legend('x','y','z')
subplot(2,1,2)
hold on
pl1 = plot(CRD_std_enu(:,1)*1000,'-b');
pl2 = plot(CRD_std_enu(:,2)*1000,'-g');
pl3 = plot(CRD_std_enu(:,3)*1000,'-r');
pl4 = plot(rmsENU(:,1)*1000,'.-b');
pl5 = plot(rmsENU(:,2)*1000,'.-g');
pl6 = plot(rmsENU(:,3)*1000,'.-r');
plot(SigmaRenu_Cov(:,1)*1000*1.8,'x-b')
plot(SigmaRenu_Cov(:,2)*1000*1.8,'x-g')
plot(SigmaRenu_Cov(:,3)*1000*1.8,'x-r')

plot([1, length(CRD_std_enu)], [rms(CRD_std_enu(:,1)), rms(CRD_std_enu(:,1))]*1000, '--b')
plot([1, length(CRD_std_enu)], [rms(CRD_std_enu(:,2)), rms(CRD_std_enu(:,2))]*1000, '--g')
plot([1, length(CRD_std_enu)], [rms(CRD_std_enu(:,3)), rms(CRD_std_enu(:,3))]*1000, '--r')
plot([1, length(CRD_std_enu)], [rms(rmsENU(:,1)), rms(rmsENU(:,1))]*1000, '--b')
plot([1, length(CRD_std_enu)], [rms(rmsENU(:,2)), rms(rmsENU(:,2))]*1000, '--g')
plot([1, length(CRD_std_enu)], [rms(rmsENU(:,3)), rms(rmsENU(:,3))]*1000, '--r')
plot(range_w, CRD_std_enu(range_w,:)*1000,'*k')
plot(range_w, rmsENU(     range_w,:)*1000,'*k')
ylabel('[mm]')
legend([pl1, pl2 pl3 pl4 pl5 pl6], {['SNX Sigma_e, rms = ' , num2str(rms(CRD_std_enu(:,1))*1000), ' mm'], ...
                                    ['SNX Sigma_n, rms = ' , num2str(rms(CRD_std_enu(:,2))*1000), ' mm'], ...
                                    ['SNX Sigma_u, rms = ' , num2str(rms(CRD_std_enu(:,3))*1000), ' mm'], ...
                                    ['PLT Sigma_e, rms = ' , num2str(rms(rmsENU(:,1))*1000), ' mm'], ...
                                    ['PLT Sigma_n, rms = ' , num2str(rms(rmsENU(:,2))*1000), ' mm'], ...
                                    ['PLT Sigma_u, rms = ' , num2str(rms(rmsENU(:,3))*1000), ' mm']} );
hold off

%%

SigmaVenu_Cov = SigmaVenu;
clc
scalePLT_SNX_std = rmsENU./CRD_std_enu;
scalePLT_SNX_sig = rmsENU./(SigmaRenu_Cov*1.8);
scaleSNX_std_sig = CRD_std_enu./(SigmaRenu_Cov*1.8);
scaleSNX_sig_RV  = SigmaRenu_Cov./SigmaVenu_Cov;

try
    close (fig2)
end
fig2 = figure(2);
hold on; grid on
plot(scalePLT_SNX_std,'-')
plot(scalePLT_SNX_sig,'.-')
plot(scaleSNX_std_sig,'x-')
plot(range_w, scalePLT_SNX_std(range_w,:), '*')
plot(range_w, scalePLT_SNX_sig(range_w,:), '*')
plot(range_w, scaleSNX_std_sig(range_w,:), '*')
legend('e','n','u')
% ylim([-100 100])
hold off  

disp('mean rms ,                  E         N         U')
disp(['CRD PLT rms           : ',num2str(mean(rmsENU,     1)*1000,'%10.4f' ), '  [mm]'])
disp(['CRD SNX sdt           : ',num2str(mean(CRD_std_enu,1)*1000,'%10.4f' ), '  [mm]'])
disp(['CRD SNX 1sigma        : ',num2str(mean(SigmaRenu_Cov,1)*1000*1.8,'%10.4f' ), '  [mm]'])
disp(['VEL SNX 1sigma        : ',num2str(mean(SigmaVenu_Cov,1)*1000*1.8,'%10.4f' ), '  [mm/yr]'])
disp('Ave. scale                 E         N         U')
disp(['PLT/SNX_std           : ', num2str(mean(scalePLT_SNX_std, 1),'%10.2f')])
disp(['PLT/SNX_sig           : ', num2str(mean(scalePLT_SNX_sig, 1),'%10.2f')])
disp(['SNX_sdt/SNX_sig       :  ',num2str(mean(scaleSNX_std_sig, 1),'%10.2f')])
disp(['CRD/VEL SNX_sig       :  ',num2str(mean(scaleSNX_sig_RV, 1),'%10.2f')])
disp(['for selected stations : ', Stations(range_w)'])
disp(['PLT/SNX_std           : ', num2str(mean(scalePLT_SNX_std(range_w,:), 1),'%10.2f')])
disp(['PLT/SNX_sig           : ', num2str(mean(scalePLT_SNX_sig(range_w,:), 1),'%10.2f')])
disp(['SNX_sdt/SNX_sig       :  ',num2str(mean(scaleSNX_std_sig(range_w,:), 1),'%10.2f')])
disp(['CRD/VEL SNX_sig       :  ',num2str(mean(scaleSNX_sig_RV(range_w,:), 1),'%10.2f')])
disp(['for selected stations : ', Stations(sel_good)'])
disp(['PLT/SNX_std           : ', num2str(mean(scalePLT_SNX_std(sel_good,:), 1),'%10.2f')])
disp(['PLT/SNX_sig           : ', num2str(mean(scalePLT_SNX_sig(sel_good,:), 1),'%10.2f')])
disp(['SNX_sdt/SNX_sig       :  ',num2str(mean(scaleSNX_std_sig(sel_good,:), 1),'%10.2f')])
disp(['CRD/VEL SNX_sig       :  ',num2str(mean(scaleSNX_sig_RV(sel_good,:), 1),'%10.2f')])
%%
Station = SiteDome_list{8};
[N, E, U, MJD, Epoch] = get_PLT_timeseries('Velocity_field/FMC_IGB_W7.PLT',Station );

plot_residuals([N, E, U], Epoch)

%% 3D plots
Lon = wrapTo180(Lon);
try
    close (fig3)
end
fig3 = figure(3);
hold on
Earth_coast(2)
xlim([-6 18]);
ylim([41 53])
% plot(Lon,Lat,'.b')
% text(Lon,Lat,Stations)
scatter(Lon,Lat,rmsENU(:,1)*1000*100)
% scatter(Lon,Lat,CRD_std_enu(:,1)*1000*1000)
% scatter(Lon,Lat,scale(:,1))








