function write_LLH(Records)
% write Residuals filed files *.txt


fileID = fopen('../MAP/LLH_CRD_all.txt', 'w');

headString = '  Long [deg],   Lat [deg],       Site Radom   \n';
formatStr = '%12.7f  %12.7f  %9s %-9s \n';

fprintf(fileID, headString);
% fprintf(headString);

disp([' ... wrinting CRD file LLH_CRD.txt'])

for i = 1:length(Records.Stations)    
    llh_crd = [Records.CRD.ENU.Est(i,[1,2])];
    fprintf(fileID, formatStr, llh_crd, Records.Stations(i,:), Records.Radoms(i,:));   
%     fprintf( formatStr, llh_crd, Records.Stations(i,:), Records.Radoms(i,:));   
end
disp('Done')
fclose(fileID);

end