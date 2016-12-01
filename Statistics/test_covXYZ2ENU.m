
clc

lat = 44.1213985;
lon = 3.5812609;

sigma_xyz = [1 1 1];

cov_xyz = [sigma_xyz(1)^2 0 0;
           0 sigma_xyz(2)^2 0;
           0 0 sigma_xyz(3)^2];


R = [          -sind(lon)             cosd(lon)       0    ;
     -sind(lat)*cosd(lon)  -sind(lat)*sind(lon)  cosd(lat) ;
      cosd(lat)*cosd(lon)   cosd(lat)*sind(lon)  sind(lat)]; 
  
cov_enu = R * cov_xyz * R';

det(R)