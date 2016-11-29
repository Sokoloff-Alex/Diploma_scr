% compute strain field from deformaion

close all
clear all
clc

%% load deformation
% DeformationField = struct2array(load('dat/Alps_deformation_0.25x0.25_no_correlaion_3.mat'));
% DeformationField = struct2array(load('dat/Alps_deformation_0.5x0.5_grid.mat'));
% DeformationField = struct2array(load('dat/Alps_deformation_1x1_grid.mat'));

load('dat/ETOPO_Alps.mat');



%%
DeformationField  = Alps_deformation;
Vel = DeformationField(:,[3,4]);
LongGrid = DeformationField(:,1);
LatGrid  = DeformationField(:,2);
Strain = getStrainMap2(DeformationField);

%%
% writeStrain2GMT(     Strain, '~/Alpen_Check/MAP/Strain/StrainField_0.25x0.25.txt')
% writeStrainSum2GMT(  Strain, '~/Alpen_Check/MAP/Strain/StrainSum.txt')
% writeStrainShear2GMT(Strain, '~/Alpen_Check/MAP/Strain/StrainShear.txt')

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
etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
% Earth_coast(2)
quiver(wrapTo180(long(Selected)),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'k', 'lineWidth',2)
quiver(wrapTo180(long(iOutliers)),lat(iOutliers),Ve_res(iOutliers)*s,Vn_res(iOutliers)*s, 0, 'lineWidth',2, 'Color',[.5 .5 .5])
text(wrapTo180(long(iOutliers)), lat(iOutliers), names(iOutliers) );
% quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'Color',[.5 .5 .5]);
plotStrainNormal(Strain, 10^7*1);
% plotStrainShear(Strain, 10^7*2);
StrainSumStack = (Strain(:,3) + Strain(:,4))*10^7;
% scatter3(Strain(:,1), Strain(:,2),  StrainSumStack)
StrainGrid = stack2grid(Strain);
StrainGridSum = squeeze(StrainGrid(:,:,3) + StrainGrid(:,:,4));
h = pcolor(StrainGrid(:,:,1), StrainGrid(:,:,2),       -StrainGridSum)
shading interp
set(h,'facealpha',.2)
xlim([-4 18])
ylim([41 53])
% colorbar
% set(0, 'DefaultFigureRenderer', 'zbuffer');
hold off



%%
clc
close all

figure('renderer','opengl')
hold on
% figure(1)
% etopo_fig = showETOPO(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo);
% Earth_coast(2)
colormap( gray(256))
hold off
hold on
% freezeColors; 

% Overlay semitransparent StrainGridSum
h = pcolor(StrainGrid(:,:,1), StrainGrid(:,:,2),       -StrainGridSum);
plotStrainNormal(Strain, 10^7*1);
colormap(jet)
shading interp
set(h,'facealpha',.5)

% overlay velocity vectors: 
quiver(wrapTo180(long(Selected)),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'k', 'lineWidth',2)
quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'Color',[.5 .5 .5]);
plotStrainNormal(Strain, 10^7*1);
xlim([-2 18])
ylim([41 53])
colorbar
hold off









