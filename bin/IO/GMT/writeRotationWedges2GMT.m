function writeRotationWedges2GMT(Strain, filename)
% function to write Shear Strain Field for GMT
%
% Alexandr Sokolov, KEG
% 28.11.2016

% GMT Format :
% -Swwedge_scale/wedge_mag.
%
%    Rotational wedges. Wedge_scale sets the size of the wedges 
%       in inches (unless c, i, or p is appended). Values are multiplied
%       by Wedge_mag before plotting. 
%       For example, setting Wedge_mag to 1.e7 works well for rotations
%       of the order of 100 nanoradians/yr. 
%       Parameters are expected to be in the following columns:
% 
%           1,2: longitude, latitude, of station (-: option interchanges order) 
%           3: rotation in radians 
%           4: rotation uncertainty in radians
%

%% save results
Wedges = [Strain(:,1:2), Strain(:,8), zeros(size(Strain(:,8)))];

fileID = fopen(filename, 'w');
fprintf(fileID, '# Rotation Wegdes, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],          Wegde[rad],           Sigma[rad] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e  \n';

for i = 1:size(Wedges,1)
   fprintf(fileID, formatStr, Wedges(i,:)); 
end
fclose(fileID);


end