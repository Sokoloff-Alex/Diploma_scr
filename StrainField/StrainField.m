% compute strain field from deformaion

close all
clear all
clc

%% load deformation
DeformationField = struct2array(load('dat/Alps_deformation_0.25x0.25_no_correlaion_3.mat'));
% DeformationField = struct2array(load('dat/Alps_deformation_0.5x0.5_grid.mat'));
% DeformationField = struct2array(load('dat/Alps_deformation_1x1_grid.mat'));

load('dat/ETOPO_Alps.mat');

Vel = DeformationField(:,[3,4]);
LongGrid = DeformationField(:,1);
LatGrid  = DeformationField(:,2);

%%
Strain = getStrainMap2(DeformationField);

%%

% writeStrain2GMT(Strain, '~/Alpen_Check/MAP/Strain/StrainField_0.25x0.25.txt')

%%
% clc
% VelocityField = [long, lat, Ve_res, Vn_res, SigmaVe*10, SigmaVn*10, zeros(size(Azim))];
% writeVelocityFieldGMT(VelocityField, names, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor.txt')

%%
clc
close all
s = 200;
figure(1)
hold on; grid on;
axis equal
% geoshow(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo, 'DisplayType', 'texturemap');
% demcmap(ETOPO_Alps.Etopo_Europe);
% cptcmap('Europe')
Earth_coast(2)
% quiver(long(Selected),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'k', 'lineWidth',2)
% quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'Color',[.5 .5 .5]);
plotStrainNormal(Strain, 10^7*1);
StrainSumStack = (Strain(:,3) + Strain(:,4))*10^7;
% scatter3(Strain(:,1), Strain(:,2),  StrainSumStack)
StrainGrid = stack2grid(Strain);
StrainGridSum = squeeze(StrainGrid(:,:,3) + StrainGrid(:,:,4));
mesh(StrainGrid(:,:,1), StrainGrid(:,:,2),       -StrainGridSum)
alpha color
% alpha scaled
% % plotStrainShear( Strain, 10^7*2); % ok
xlim([-2 18])
ylim([41 53])
colorbar
set(0, 'DefaultFigureRenderer', 'zbuffer');
hold off

%%
clc
close all

figure('renderer','opengl')
hold on
% figure(1)
% etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
Earth_coast(2)
colormap( gray(256))
hold off
hold on
% freezeColors; 

% Overlay semitransparent ice speed:
h = pcolor(StrainGrid(:,:,1), StrainGrid(:,:,2),       -StrainGridSum);
colormap( jet)
shading interp
set(h,'facealpha',.2)

% overlay velocity vectors: 
quiver(long(Selected),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'k', 'lineWidth',2)
quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'Color',[.5 .5 .5]);
plotStrainNormal(Strain, 10^7*1);
xlim([-2 18])
ylim([41 53])
colorbar
hold off



%% try for GMT

StrainFieldSum = [Strain(:,1), Strain(:,2), StrainSumStack];

filename = '~/Alpen_Check/MAP/Strain/StrainSum.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, '# Normal Strain Field, Sum of E1 + E2, \n');
fprintf(fileID, '# Long [deg],   Lat [deg],   Eps1+Eps2 [strain/yr] \n');
formatStr = '%12.7f  %12.7f  %+12e  \n';

for i = 1:size(Strain,1)
   fprintf(fileID, formatStr, StrainFieldSum(i,:)); 
end
fclose(fileID);







