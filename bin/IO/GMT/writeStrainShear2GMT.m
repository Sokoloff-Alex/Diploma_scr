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

Alpha = Strain(:,7);

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

%% save in 2 files, depending on rotation wedge w : CW & CCW
% to plot with -Sr
% long, lat , Ve, Vn, 0,0,0
% prepare CW and CCW

w = Strain(:,8);
N = 1:length(w);
nCCW = length(w(w >= 0));
nCW  = length(w) - nCCW;

ShearStrainCCW = [Strain(N(w >= 0),1:2), Strain(N(w >= 0),6)*10^9, zeros(nCCW,1), -Strain(N(w >= 0),7)  ];
ShearStrainCW  = [Strain(N(w <  0),1:2), Strain(N(w <  0),6)*10^9, zeros(nCW, 1), -Strain(N(w <  0),7)  ];

ShearArrowCCW = zeros(nCCW,4);
ShearArrowCW  = zeros(nCW,4);

for i = 1:nCCW
   ShearArrowCCW(i,:) = [ShearStrainCCW(i,1:2), ShearStrainCCW(i,3)*cosd(-ShearStrainCCW(i,5)), ShearStrainCCW(i,3)*sind(-ShearStrainCCW(i,5)) ]; 
end

for i = 1:nCW
   ShearArrowCW(i,:)  = [ShearStrainCW(i,1:2),  ShearStrainCW(i,3)*cosd(-ShearStrainCW(i,5)),   ShearStrainCW(i,3)*sind(-ShearStrainCW(i,5)) ]; 
end

%%
% write CCW, +w
fileID = fopen([filename(1:end-4),'_CCW.txt'], 'w');
fprintf(fileID, '# Shear Strain Field, w>=0  => CCW, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],         Eps1_e [nstrain/yr],    Eps1_n [nstrain/yr] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e  \n';
for i = 1:nCCW
   fprintf(fileID, formatStr,  ShearArrowCCW(i,:)); 
   fprintf(fileID, formatStr, [ShearArrowCCW(i,1:2), -ShearArrowCCW(i,3:4) ]); 
end
fclose(fileID);

% write CW, -w
fileID = fopen([filename(1:end-4),'_CW.txt'], 'w');
fprintf(fileID, '# Shear Strain Field, w>0  => CW, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],         Eps1_e [nstrain/yr],    Eps1_n [nstrain/yr] \n');
formatStr = '%12.7f  %12.7f  %+20e  %+20e \n';
for i = 1:nCW
   fprintf(fileID, formatStr, ShearArrowCW(i,:)); 
   fprintf(fileID, formatStr, [ShearArrowCW(i,1:2), -ShearArrowCW(i,3:4) ]); 
end
fclose(fileID);



end