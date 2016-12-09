function writeVelocityFieldVertical2GMT(Long, Lat, Vu, filename)
% function to write Velocity Field Vertical for GMT
%
% Alexandr Sokolov, KEG
% 09.12.2016

VelocityField = [Long, Lat, Vu];

% filename = '~/Alpen_Check/MAP/VelocityField/VelocityFieldVertical.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '# Velocity Field Vertical, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],      Vu [mm/yr] \n');
formatStr = '%12.7f  %12.7f  %12.5f  \n';

for i = 1:size(VelocityField,1)
   fprintf(fileID, formatStr, VelocityField(i,:)); 
end
fclose(fileID);


end