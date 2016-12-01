function [CovENU, SigmaENU] = covXYZ2ENU(CovXYZ, lat ,lon)
% function to transform Covariance matrix from xyz to 
% local geodetical ENU system
%
% According to https://www.ngs.noaa.gov/CORS/Articles/SolerChin1985.pdf
% 
% Cenu = R  * Cxyz * R'
% Cxyz = R' * Cenu * R 
%
% Alexandr Sokolov, KEG
% 01.12.2016

% rotation marix

if (size(CovXYZ,1) == 1) || (size(CovXYZ,2) == 1) 
   SigmaXYZ =  CovXYZ;
   
   CovXYZ = [SigmaXYZ(1)^2    0               0       ;
              0           SigmaXYZ(2)^2       0       ;
              0               0         SigmaXYZ(3)^2];
end

R = [          -sind(lon)             cosd(lon)       0    ;
     -sind(lat)*cosd(lon)  -sind(lat)*sind(lon)  cosd(lat) ;
      cosd(lat)*cosd(lon)   cosd(lat)*sind(lon)  sind(lat)];

CovENU = R * CovXYZ * R';
  
SigmaENU = sqrt([CovENU(1,1), CovENU(2,2), CovENU(3,3)]);


   

  
end