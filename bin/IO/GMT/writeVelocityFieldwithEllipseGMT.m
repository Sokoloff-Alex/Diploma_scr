function writeVelocityFieldwithEllipseGMT(VelocityField, SiteNames, filename)
% function to wtrite Velocity firld for GMT
% 
% Alexandr Sokolov, KEG
% 30.11.2016

% format :
% -Srvelscale/confidence/fontsize
% 
%     Velocity ellipses in rotated convention. 
%         Parameters are expected to be in the following columns:
% 
%             1,2: longitude, latitude, of station (-: option interchanges order) 
%             3,4: eastward, northward velocity (-: option interchanges order) 
%             5,6: semi-major, semi-minor axes 
%             7: COUNTER-CLOCKWISE angle, in degrees, from HORIZONTAL axis to MAJOR axis of ellipse. 
%             8: name of station (optional)

fileID = fopen(filename, 'w');
fprintf(fileID, '# Velocity Field, \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma A [m/yr], Sigma B [m/yr],      Angle [deg],   Site \n');

formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f %9s \n';

for i = 1:size(VelocityField,1)
   data = [VelocityField(i,1:4), max(VelocityField(i,5:6)), min(VelocityField(i,5:6)), VelocityField(i,7) ]; 
   fprintf(fileID, formatStr, data, SiteNames{i}); 
end
fclose(fileID);

end