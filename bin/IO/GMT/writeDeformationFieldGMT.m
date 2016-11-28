function writeDeformationFieldGMT(DeformationField, filename)


fileID = fopen(filename, 'w');

fprintf(fileID, '# Velocity Field, \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr],      Angle [deg] \n');


formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f\n';

for i = 1:size(DeformationField,1)
   fprintf(fileID, formatStr, DeformationField(i,:)); 
end
fclose(fileID);

end