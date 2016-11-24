function write_velocity_field_files(Records)
% write velocity filed files *.txt


fileID_h = fopen('../Velocity_field_horizontal.txt', 'w');
fileID_v = fopen('../Velocity_field_vertical.txt',   'w');
fileID_u = fopen('../Velocity_field_vertical_Uplift.txt',   'w');
fileID_s = fopen('../Velocity_field_vertical_subduction.txt',   'w');

fprintf(fileID_h, '# Velocity field, horizontal component \n');
fprintf(fileID_v, '# Velocity field, vertical component \n');
fprintf(fileID_u, '# Velocity field, vertical component, uplift only \n');
fprintf(fileID_s, '# Velocity field, vertical component, subsidence only \n');


headString = '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr], Angle [deg],   Site \n';
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


end