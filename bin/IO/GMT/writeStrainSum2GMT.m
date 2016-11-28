function writeStrainSum2GMT(Strain, filename)
% function to write Shear Strain Field for GMT
%
% Alexandr Sokolov, KEG
% 28.11.2016

StrainSumStack = (Strain(:,3) + Strain(:,4))*10^9; % [nstrain/yr]
StrainFieldSum = [Strain(:,1), Strain(:,2), StrainSumStack];

% filename = '~/Alpen_Check/MAP/Strain/StrainSum.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '# Normal Strain Field, Sum of (E1 + E2), \n');
fprintf(fileID, '# Long [deg],   Lat [deg],   Eps1+Eps2 [nstrain/yr] \n');
formatStr = '%12.7f  %12.7f  %+12e  \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, StrainFieldSum(i,:)); 
end
fclose(fileID);


end