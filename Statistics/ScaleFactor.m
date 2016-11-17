%%  Determine Scale factor 

close all
clear all
clc

%%
SINEX = readSNX('STA/FMC_IGB_W7.SNX');
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
[rmsENU(:,2), rmsENU(:,1), rmsENU(:,3), Sites] = get_PLT_residuals('STA/FMC_IGB_W7.PLT', SiteDome_list);

%%
find(strcmp(SiteDome_list,'GRAS GRASigb03'))
range_w = [101 102 103 149 160 291 292 296];

%%
[CRD_std_enu(:,1), CRD_std_enu(:,2), CRD_std_enu(:,3)] = XYZ2ENU(CRD, CRD_std);
try
    close (fig1)
end
fig1 = figure(1);
subplot(2,1,1)
hold on
plot(CRD_std*1000)
ylabel('[mm]')
legend('x','y','z')
subplot(2,1,2)
hold on
plot( CRD_std_enu(:,1)*1000,'-b')
plot(-CRD_std_enu(:,2)*1000,'-g')
plot( CRD_std_enu(:,3)*1000,'-r')
plot(rmsENU(:,1)*1000,'.-b')
plot(rmsENU(:,2)*1000,'.-g')
plot(rmsENU(:,3)*1000,'.-r')
plot(range_w, CRD_std_enu(range_w,:)*1000,'*k')
plot(range_w, rmsENU(     range_w,:)*1000,'*k')
ylabel('[mm]')
legend('e','n','u')
hold off

%%
clc
scale = abs(rmsENU ./ [CRD_std_enu(:,1) -CRD_std_enu(:,2), CRD_std_enu(:,3) ]);

try
    close (fig2)
end
fig2 = figure(2);
hold on; grid on
plot(scale)
plot(range_w, scale(range_w,:), '*')
legend('e','n','u')
% ylim([-100 100])
hold off  
disp('mean rms       E         N         U')
disp(['SNX        : ',num2str(mean(CRD_std_enu,1)*1000,'%10.3f' )])
disp(['PLT        : ',num2str(mean(rmsENU,     1)*1000,'%10.3f' )])
disp('Ave. scale       E         N         U')
disp(['PLT/SNX all: ', num2str(mean(scale,           1),'%10.2f')])
disp(['PLT/SNX sel: ', num2str(mean(scale(range_w,:),1),'%10.2f')])

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








