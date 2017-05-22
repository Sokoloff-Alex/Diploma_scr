function writeVelocityFieldVertical(Long, Lat, Vu, Site, filename)
% function to write Velocity Field Vertical for GMT
%
% Alexandr Sokolov, KEG
% 09.12.2016

VelocityField = [Long, Lat, zeros(size(Vu)),  Vu];

% filename = '~/Alpen_Check/MAP/VelocityField/VelocityFieldVertical_out.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '# Velocity Field Vertical, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],       Ve(=0)       Vu [mm/yr]  Site \n');
formatStr = '%12.7f  %12.7f  %12.5f %12.5f     %s  \n';

for i = 1:size(VelocityField,1)
   fprintf(fileID, formatStr, VelocityField(i,:), cell2mat(Site(i))); 
end
fclose(fileID);

end