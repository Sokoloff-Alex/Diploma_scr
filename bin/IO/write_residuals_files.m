function write_residuals_files(Records)
% write Residuals filed files *.txt


fileID_h = fopen('../MAP/Residuals_field_horizontal.txt', 'w');
fileID_v = fopen('../MAP/Residuals_field_vertical.txt',   'w');
fileID_3d = fopen('../MAP/Residuals_field_3D.txt',   'w');

headString = '  Long [deg],   Lat [deg],     RMS E [mm],  RMS N [mm],   Site \n';
formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %9s \n';

fprintf(fileID_h, headString);
fprintf(fileID_v, headString);

disp([' ... Residuals vector field files *.txt'])
fileID_i = fopen('../Residuals_field_info.txt', 'w');
fprintf(fileID_i, ['4 52.20         15 0 0 5  Residuals. ']);


for i = 1:length(Records)
    
    data_h = [Records.CRD.ENU.Est(i,[1,2]),     Records.CRD.ENU.rms(i,[1,2])*1000];
    data_v = [Records.CRD.ENU.Est(i,[1,2]), 0,  Records.CRD.ENU.rms(i,3)*1000];
    
    rms_3d = sqrt(Records.CRD.ENU.rms(i,:)*Records.CRD.ENU.rms(i,:)')*1000;
    data_rms_3D = [Records.CRD.ENU.Est(i,[1,2]), 0, rms_3d ];
 
 
    fprintf(fileID_h, formatStr, data_h, Records.Stations(i,:));
    fprintf(fileID_v, formatStr, data_v, Records.Stations(i,:));
    fprintf(fileID_3d, formatStr, data_rms_3D, Records.Stations(i,:));
    
    
end
disp('Done')
fclose(fileID_h);
fclose(fileID_v);
fclose(fileID_i);
fclose(fileID_3d);


end