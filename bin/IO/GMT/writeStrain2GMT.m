function writeStrain2GMT(Strain, filename)
% function to write Strain Field for GMT
%
% Alexandr Sokolov, KEG
% 24.11.2016

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

%% save results

scale = 10^9; % [strain/yr] -> [nanostrain/yr]

StrainField = [Strain(:,1:2), Strain(:,3:4)*scale, -Strain(:,5) ];

fileID = fopen(filename, 'w');
fprintf(fileID, '# Normal Strain Field, in [nstrain/yr] \n');
fprintf(fileID, '# Long [deg],   Lat [deg],        Eps1 [nstrain/yr],    Eps2 [nstrain/yr],  AziEps2 [deg] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e   %10.2f \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, StrainField(i,:)); 
end
fclose(fileID);

%% split into Extention and Compression

StrainFieldDilataion = zeros(size(StrainField));
StrainFieldCompression = zeros(size(StrainField));

for i = 1:size(StrainField,1)
    L1 = StrainField(i,3);
    L2 = StrainField(i,4);
    if L1 < 0; L1 = 0; end
    if L2 < 0; L2 = 0; end
    StrainFieldDilataion(i,:)  = [StrainField(i,1:2),L1, L2, StrainField(i,5) ];
end

for i = 1:size(StrainField,1)
    L1 = StrainField(i,3);
    L2 = StrainField(i,4);
    if L1 > 0; L1 = 0; end
    if L2 > 0; L2 = 0; end
    StrainFieldCompression(i,:) = [StrainField(i,1:2),L1, L2, StrainField(i,5)];
end

%% write Dilatation
fileID = fopen([filename(:,1:end-4),'_Dilatation.txt'], 'w');
fprintf(fileID, '# Normal Strain Field, Dilatation only \n');
fprintf(fileID, '# Long [deg],   Lat [deg],        Eps1 [nstrain/yr],    Eps2 [nstrain/yr],  AziEps2 [deg] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e   %10.2f \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, StrainFieldDilataion(i,:)); 
end
fclose(fileID);

%% Write Compression
fileID = fopen([filename(:,1:end-4),'_Compression.txt'], 'w');
fprintf(fileID, '# Normal Strain Field, Compression only \n');
fprintf(fileID, '# Long [deg],   Lat [deg],        Eps1 [nstrain/yr],    Eps2 [nstrain/yr],  AziEps2 [deg] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e   %10.2f \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, StrainFieldCompression(i,:)); 
end
fclose(fileID);


end