
clear all
% close all
clc

cd ~/Alpen_Check/Matlab_dir/
%%

[Stations, Radoms, Records] = parse_OUT_file('/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FCS_out4.OUT');
% [Stations, Radoms, Records] = parse_OUT_file('Results/FCSEPN.OUT');

%% write velocity filed files *.txt
close all

fileID_h = fopen('../Velocity_field_horizontal.txt', 'w');
fileID_v = fopen('../Velocity_field_vertical.txt',   'w');
fileID_u = fopen('../Velocity_field_vertical_Uplift.txt',   'w');
fileID_s = fopen('../Velocity_field_vertical_subduction.txt',   'w');

headString = '  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr], Angle [deg],   Site \n';
formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %12.5f  %12.5f  %16.5f %9s \n';

fprintf(fileID_h, headString);
fprintf(fileID_v, headString);
fprintf(fileID_u, headString);
fprintf(fileID_s, headString);

Frame = 'ETRF';
disp([' ... writing vector field files *.txt in ', Frame, ' Frame'])
fileID_i = fopen('../Velocity_field_info.txt', 'w');
fprintf(fileID_i, ['4 52.20         15 0 0 5  ', Frame, ' . ']);


for i = 1:Records.NumberOfStations 
    if strcmp(Frame,'ETRF')
        [R_E, V_E]   = ITRF2ETRF(Records.CRD.XYZ.Est(i,:),Records.VEL.XYZ.Est(i,:));
        [Ve, Vn, Vu, lat, lon, h] = XYZ2ENU(R_E, V_E);
        data_h = [lon, lat, Ve, Vn, Records.CRD.ENU.Ellipse.Sigmas(i,:),       Records.CRD.ENU.Ellipse.Angle(i)];
        data_v = [lon, lat, 0,  Vu, Records.CRD.ENU.Ellipsoid.Sigmas(i,[3,3]), 0];
    elseif strcmp(Frame, 'ITRF')
        data_h = [Records.CRD.ENU.Est(i,[1,2]), Records.VEL.ENU.Est(i,[1,2]), Records.CRD.ENU.Ellipse.Sigmas(i,:),   Records.CRD.ENU.Ellipse.Angle(i)];
        data_v = [Records.CRD.ENU.Est(i,[1,2]), 0, Records.VEL.ENU.Est(i,3),  Records.CRD.ENU.Ellipsoid.Sigmas(i,[3,3]), 0];
    else
        disp('choose b/w: ETRF or ITRF')
        break
    end
 
    fprintf(fileID_h, formatStr, data_h, Records.Stations(i,:));
    fprintf(fileID_v, formatStr, data_v, Records.Stations(i,:));
    if data_v(4) >= 0
        fprintf(fileID_u, formatStr, data_v, Records.Stations(i,:));
    else
        fprintf(fileID_s, formatStr, data_v, Records.Stations(i,:));
    end
    
end
disp('Done')
fclose(fileID_h);
fclose(fileID_v);
fclose(fileID_u);
fclose(fileID_s);
fclose(fileID_i);

%%
close all
scale = 1000*1000*10;
clr = lines(5);
figure(1)
hold on
grid on
axis equal
title('Velocity field in ITRF')
pl1 = quiver3(Records.CRD.XYZ.Apr(:,1),Records.CRD.XYZ.Apr(:,2),Records.CRD.XYZ.Apr(:,3),  Records.VEL.XYZ.Apr(:,1)*scale, Records.VEL.XYZ.Apr(:,2)*scale, Records.VEL.XYZ.Apr(:,3)*scale ,0,'Color',clr(1,:),'LineWidth', 2)
pl2 = quiver3(Records.CRD.XYZ.Est(:,1),Records.CRD.XYZ.Est(:,2),Records.CRD.XYZ.Est(:,3),  Records.VEL.XYZ.Est(:,1)*scale, Records.VEL.XYZ.Est(:,2)*scale, Records.VEL.XYZ.Est(:,3)*scale ,0,'Color',clr(2,:),'LineWidth', 2)
%     quiver3(Records.CRD.XYZ.Est(:,1),Records.CRD.XYZ.Est(:,2),Records.CRD.XYZ.Est(:,3),  Records.VEL.XYZ.Cor(:,1)*scale, Records.VEL.XYZ.Cor(:,2)*scale, Records.VEL.XYZ.Cor(:,3)*scale ,0,'Color',clr(3,:),'LineWidth', 2)
CRD3 = Records.CRD.XYZ.Apr + Records.VEL.XYZ.Apr*scale;
pl3 = quiver3(CRD3(:,1),CRD3(:,2),CRD3(:,3),  Records.VEL.XYZ.Cor(:,1)*scale, Records.VEL.XYZ.Cor(:,2)*scale, Records.VEL.XYZ.Cor(:,3)*scale ,0,'Color',clr(3,:),'LineWidth', 2)
legend([pl1, pl2, pl3], 'A priory', 'Estimate', 'Correction')

%% in ETRF 

[CRD_Apr_ETRF, VEL_Apr_ETFR]   = ITRF2ETRF(Records.CRD.XYZ.Apr, Records.VEL.XYZ.Apr);
[CRD_Est_ETRF, VEL_Est_ETFR]   = ITRF2ETRF(Records.CRD.XYZ.Est, Records.VEL.XYZ.Est);

CRD3 = Records.CRD.XYZ.Apr + Records.VEL.XYZ.Apr*scale;
%%

close all
scale = 1000*1000*50;
clr = lines(5);
figure(1)
hold on
grid on
axis equal
title('Velocity field in ETRF')
pl1 = quiver3(CRD_Apr_ETRF(:,1),CRD_Apr_ETRF(:,2),CRD_Apr_ETRF(:,3),  VEL_Apr_ETFR(:,1)*scale, VEL_Apr_ETFR(:,2)*scale, VEL_Apr_ETFR(:,3)*scale ,0,'Color',clr(1,:),'LineWidth', 2)
pl2 = quiver3(CRD_Est_ETRF(:,1),CRD_Est_ETRF(:,2),CRD_Est_ETRF(:,3),  VEL_Est_ETFR(:,1)*scale, VEL_Est_ETFR(:,2)*scale, VEL_Est_ETFR(:,3)*scale ,0,'Color',clr(2,:),'LineWidth', 2)
% quiver3(Records.CRD.XYZ.Est(:,1),Records.CRD.XYZ.Est(:,2),Records.CRD.XYZ.Est(:,3),  Records.VEL.XYZ.Cor(:,1)*scale, Records.VEL.XYZ.Cor(:,2)*scale, Records.VEL.XYZ.Cor(:,3)*scale ,0,'Color',clr(3,:),'LineWidth', 2)
CRD3 = CRD_Apr_ETRF + VEL_Apr_ETFR*scale;
pl3 = quiver3(CRD3(:,1),CRD3(:,2),CRD3(:,3),  Records.VEL.XYZ.Cor(:,1)*scale, Records.VEL.XYZ.Cor(:,2)*scale, Records.VEL.XYZ.Cor(:,3)*scale ,0,'Color',clr(3,:),'LineWidth', 2)
legend([pl1, pl2, pl3], 'A priory', 'Estimate', 'Correction')



%%

scale = 1000;
clr = lines(5);

fig1 = figure(4);
hold on
grid on
% grid minor
% Earth_coast(2)
title('Velocity field, vertical component')
text(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2), Stations);
plot(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2), '.')
% quiver(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2),  Records.VEL.ENU.Apr(:,1)*scale, Records.VEL.ENU.Apr(:,2)*scale, 0,'Color',clr(1,:),'LineWidth', 2)
% plt1 = quiver(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2),  zeros(size(Records.VEL.ENU.Apr(:,1))), Records.VEL.ENU.Apr(:,3)*scale, 0,'Color',clr(1,:),'LineWidth', 2);
% plt2 = quiver(Records.CRD.ENU.Est(:,1),Records.CRD.ENU.Est(:,2),  zeros(size(Records.VEL.ENU.Est(:,1))), Records.VEL.ENU.Est(:,3)*scale, 0,'LineWidth', 2);
% plt3 = quiver(Records.CRD.ENU.Est(:,1),Records.CRD.ENU.Est(:,2),  zeros(size(Records.VEL.ENU.Cor(:,1))), Records.VEL.ENU.Cor(:,3)*scale, 0,'Color',clr(3,:),'LineWidth', 2);

xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
hold on
grid on
% legend([plt1, plt2, plt3], 'A priory', 'Estimated', 'Correction')

for i = 1:Records.NumberOfStations
    plot(Records.CRD.ENU.Apr(i,1),Records.CRD.ENU.Apr(i,2),'.b')
%     ellipce_2D(Records.CRD.ENU.Ellipse.Sigmas(i,:),       Records.CRD.ENU.Ellipse.Angle(i), Records.CRD.ENU.Est(i,1:2), scale); % horzontal
    ellipce_2D(Records.CRD.ENU.Ellipsoid.Sigmas(i,[3,3]), 0,                                Records.CRD.ENU.Est(i,1:2), scale); % vertical
    if (Records.VEL.ENU.Est(i,3) >= 0)
        plt2 = quiver(Records.CRD.ENU.Est(i,1),Records.CRD.ENU.Est(i,2),  0, Records.VEL.ENU.Est(i,3)*scale, 0,'Color',[0 0.5 0], 'LineWidth', 2);
    else
        plt2 = quiver(Records.CRD.ENU.Est(i,1),Records.CRD.ENU.Est(i,2),  0, Records.VEL.ENU.Est(i,3)*scale, 0,'Color','r',  'LineWidth', 2);
    end
end

% plt2 = quiver(Records.CRD.ENU.Est(:,1),Records.CRD.ENU.Est(:,2),  zeros(size(Records.VEL.ENU.Est(:,1))), Records.VEL.ENU.Est(:,3)*scale,0,'LineWidth', 2);
xlim([-5 17])
ylim([41 53])


%%
close all
figure(3)
hold on
grid on
title('Velocity field, horizontal component')

plot(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2), '.b')
% quiver(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2),  Records.VEL.ENU.Apr(:,1)*scale, Records.VEL.ENU.Apr(:,2)*scale, 0,'Color',clr(1,:),'LineWidth', 2)
% plt1 = quiver(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2),  Records.VEL.ENU.Apr(:,1)*scale, Records.VEL.ENU.Apr(:,2)*scale, 0,'Color',clr(1,:),'LineWidth', 2)
% plt2 = quiver(Records.CRD.ENU.Est(:,1),Records.CRD.ENU.Est(:,2),  Records.VEL.ENU.Est(:,1)*scale, Records.VEL.ENU.Est(:,2)*scale, 0,'Color',clr(2,:),'LineWidth', 2)
plt3 = quiver(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2),  Records.VEL.ENU.Cor(:,1)*scale, Records.VEL.ENU.Cor(:,2)*scale, 0,'Color',clr(3,:),'LineWidth', 2)
% text(Records.CRD.ENU.Apr(:,1),Records.CRD.ENU.Apr(:,2), Stations);
xlabel('Longitude, [deg]')
ylabel('Latitude,  [deg]')
legend([plt1, plt2, plt3], 'A priory', 'Estimated', 'Correction')

