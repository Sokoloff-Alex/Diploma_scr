function writeDeformationFieldGMT(DeformationField, filename , flag)


fileID = fopen(filename, 'w');

fprintf(fileID, '# Deformation Field, \n');
if ismember(flag, {'-e', 'Ellipse', 'azim', 'angle'})
    fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr],      Angle [deg] \n');
    formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f\n';
elseif ismember(flag, {'Corr', '-c'})
    fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr],      CorrEN      \n');
    formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f\n';
elseif ismember(flag, {'Simple', '-s','noCov'})
    fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr] \n');
    formatStr = '%12.7f  %12.7f  %12.5f  %12.5f  \n';
end

for i = 1:size(DeformationField,1)
   fprintf(fileID, formatStr, DeformationField(i,:)); 
end
fclose(fileID);

disp('Done')

end