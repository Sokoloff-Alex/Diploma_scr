% test Cov_xyz2enu
%  AIGL 10059M003        X      4578299.94982    4578299.95053       0.00071       0.00013
%                        Y       286538.94753     286538.94759       0.00007       0.00004
%                        Z      4418911.81466    4418911.81550       0.00084       0.00012
% 
%                        U         1618.77252       1618.77362       0.00110       0.00016     0.00016    2.4
%                        N         44.1213985       44.1213985       0.00011       0.00005     0.00004   96.8     0.00004   96.8
%                        E          3.5812609        3.5812609       0.00002       0.00004     0.00005    0.0     0.00005
% 
%  AIGL 10059M003        VX          -0.01208         -0.01196       0.00012       0.00002
%                        VY           0.01876          0.01884       0.00008       0.00000
%                        VZ           0.01205          0.01209       0.00004       0.00001
% 
%                        VU           0.00058          0.00069       0.00012       0.00002     0.00002    2.8
%                        VN           0.01623          0.01617      -0.00005       0.00001     0.00000   97.4     0.00000   97.2
%                        VE           0.01948          0.01955       0.00007       0.00000     0.00001    0.1     0.00001

clc
format short
lat1 = 44.1213985;
long1 = 3.5812609;
sigma_xyz = [0.00013  0.00004  0.00012]*1000;

sigma_xyz = [1 2 3];

cov_xyz = [sigma_xyz(1)^2 0 0;
           0 sigma_xyz(2)^2 0;
           0 0 sigma_xyz(3)^2];
   
[cov_enu, sigma_enu] = covXYZ2ENU(cov_xyz, lat1, long1 );
corr_en = cov_enu(1,2) / (sigma_enu(1) * sigma_enu(2));
angle = 90 + 1/2 * atand(2*cov_enu(1,2) / (cov_enu(1,1) - cov_enu(2,2)));
if sigma_enu(1) < sigma_enu(1);
    angle = angle + 90;
end
azim = 90 - angle;

% reverse
cov_xyz2 = covENU2XYZ(cov_enu, lat1, long1);

[cov_xyz cov_enu cov_xyz2]
[det(cov_xyz), det(cov_enu), det(cov_xyz2)]

%%
close all
figure(1)
hold on
grid on
axis equal
error_ellipse([cov_enu(1:2, 1:2)], [0 0 ],0.67, 'r')
% error_ellipse([cov_xyz(1:2, 1:2)], [0 0 ],0.67, 'b')
