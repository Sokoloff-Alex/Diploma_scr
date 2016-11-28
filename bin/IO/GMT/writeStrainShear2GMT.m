function writeStrainShear2GMT(Strain, filename)
% function to write Shear Strain Field for GMT
%
% Alexandr Sokolov, KEG
% 28.11.2016

% GMT Format :
% -Sxcross_scale
%       gives Strain crosses. Cross_scale sets the size 
%       of the cross in inches (unless c, i, or p is appended). 
%       Parameters are expected to be in the following columns:
% 
%         1,2: longitude, latitude, of station (-: option interchanges order) 
%         3: eps1, the most extensional eigenvalue of strain tensor, 
%           with extension taken positive. 
%         4: eps2, the most compressional eigenvalue of strain tensor, 
%           with extension taken positive. 5: azimuth of eps2
%           in degrees CW from North.
%
% output only Eps1 = +(sqrt( 1/4*(exx-eyy)^2 + exy^2 )) Max value,
% Eps2 = -Eps1, not interesting, in outout file Eps2 = 0
%

%% save results
ShearStrainField = [Strain(:,1:2), Strain(:,6)*10^9, zeros(size(Strain(:,6))), -Strain(:,7) ];

fileID = fopen(filename, 'w');
fprintf(fileID, '# Shear Strain Field, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],         Eps1 [nstrain/yr],    Eps2 =-Eps1 (=0!), AziEps2 [deg] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e   %10.2f \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, ShearStrainField(i,:)); 
end
fclose(fileID);


end