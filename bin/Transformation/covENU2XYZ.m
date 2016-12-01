function [CovXYZ, SigmaXYZ] = covENU2XYZ(CovENU, lat ,lon)
% function to transform Covariance matrix from local geodetical ENU system
% to xyz
% 
%
% According to https://www.ngs.noaa.gov/CORS/Articles/SolerChin1985.pdf
% 
% Cenu = R  * Cxyz * R'
% Cxyz = R' * Cenu * R 
%
% Alexandr Sokolov, KEG
% 01.12.2016

% rotation marix

if (size(CovENU,1) == 1) || (size(CovENU,2) == 1) 
   SigmaENU =  CovENU;
   
   CovENU = [SigmaENU(1)^2    0               0       ;
              0           SigmaENU(2)^2       0       ;
              0               0         SigmaENU(3)^2];
end

R = [          -sind(lon)             cosd(lon)       0    ;
     -sind(lat)*cosd(lon)  -sind(lat)*sind(lon)  cosd(lat) ;
      cosd(lat)*cosd(lon)   cosd(lat)*sind(lon)  sind(lat)]'; % <<< TRANSFORMED!

CovXYZ = R * CovENU * R';
  
SigmaXYZ = sqrt([CovXYZ(1,1), CovXYZ(2,2), CovXYZ(3,3)]);

  
end