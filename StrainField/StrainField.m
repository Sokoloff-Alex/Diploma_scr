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

writeStrain2GMT(Strain, '~/Alpen_Check/MAP/Strain/StrainField_0.25x0.25.txt')

%%
clc
VelocityField = [long, lat, Ve_res, Vn_res, SigmaVe*10, SigmaVn*10, zeros(size(Azim))];


writeVelocityFieldGMT(VelocityField, names, '~/Alpen_Check/MAP/VelocityField/VelocityField_hor.txt')



%%




%%

clc
%%
close all
s = 200;
figure(1)
hold on; grid on;
axis equal
geoshow(ETOPO_Alps.Etopo_Europe, ETOPO_Alps.refvec_Etopo, 'DisplayType', 'texturemap');
demcmap(ETOPO_Alps.Etopo_Europe);
cptcmap('Europe')
% quiver(long(Selected),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'k', 'lineWidth',2)
% quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'Color',[.5 .5 .5]);
plotStrainNormal(Strain, 10^7*1);
% % plotStrainShear( Strain, 10^7*2); % ok
xlim([-2 18])
ylim([41 53])
set(0, 'DefaultFigureRenderer', 'zbuffer');
hold off






