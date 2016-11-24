function writeVelocityFieldGMT(VelocityField, SiteNames, filename)


fileID = fopen(filename, 'w');

fprintf(fileID, '# Velocity Field, \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma E [m/yr], Sigma N [m/yr],      Angle [deg],   Site \n');


formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f %9s \n';

for i = 1:size(VelocityField,1)
   fprintf(fileID, formatStr, VelocityField(i,:), SiteNames{i}); 
end
fclose(fileID);

end