function writeVelocityFieldwithCovGMT(VelocityField, SiteNames, filename)
% function to wtrite Velocity firld for GMT
% 
% Alexandr Sokolov, KEG
% 30.11.2016

% format :
% -Sevelscale/confidence/fontsize.
% 
%     Velocity ellipses in (N,E) convention. 
%       Parameters are expected to be in the following columns:
%         1,2: longitude, latitude of station (-: option interchanges order) 
%         3,4: eastward, northward velocity (-: option interchanges order) 
%         5,6: uncertainty of eastward, northward velocities (1-sigma) (-: option interchanges order) 
%         7: correlation between eastward and northward components 
%         8: name of station (optional).



fileID = fopen(filename, 'w');
fprintf(fileID, '# Velocity Field, \n');
fprintf(fileID, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma Ve [m/yr], Sigma Vn [m/yr],      CorrVen,     Site \n');

fileID2 = fopen([filename(1:end-4),'_TH.txt'], 'w');
fprintf(fileID2, '# Velocity Field, \n');
fprintf(fileID2, '#  Long [deg],   Lat [deg],     Vel E [m/yr],  Vel N [m/yr], Sigma Ve [m/yr], Sigma Vn [m/yr],      CorrVen,     Site \n');


formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %15e  %15e  %16.5f %9s \n';

for i = 1:size(VelocityField,1)
    SigE = (VelocityField(i,5));
    SigN = (VelocityField(i,6));
    
    SigmaTH = 0.0001;
    
    if SigE < SigmaTH || SigN < SigmaTH
       if SigE < SigmaTH
          SigE = SigmaTH;
       end
       if SigN < SigmaTH
          SigN = SigmaTH;
       end
       data = [VelocityField(i,1:4), SigE, SigN, VelocityField(i,7)]; 
       fprintf(fileID2, formatStr, data, SiteNames{i}); 
%        fprintf(fileID, formatStr, data, SiteNames{i}); 
    else
       data = [VelocityField(i,1:4), SigE, SigN, VelocityField(i,7)]; 
       fprintf(fileID, formatStr, data, SiteNames{i}); 
    end
end
fclose(fileID);
fclose(fileID2);

end