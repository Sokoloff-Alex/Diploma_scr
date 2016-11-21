% script to validate SCL result with other velocity fields

%% make structure from France Supl.Data


close all
clear all
clc

%%
clc
% data = Francevelocityfieldsupldata;
% 
% FieldNames = {'Sites','Lat','Long','Ve','Vn','Vu', ...     
%               'SigmaVe','SigmaVn','SigmaVu', ...
%               'Period','Days','years', ...
%               'kappa_E','kappa_N','kappa_U',};
% 
% FranceVELStruct = cell2struct(data(1:end,:)', FieldNames);
% FranceVELtable = struct2table(FranceVELStruct);
% save('FranceVELtable.mat', 'FranceVELtable');

load('Validation/FranceVELtable.mat')
FRvel = FranceVELtable;
clear FranceVELtable;

%% Interpolate ALP_NET velocities at France sites, to compare
clc
tic
p = length(FRvel.Lat);
range = 1:length(lat);
Max_Dist = 150; % km
lim = 10;
for i = 1:p
iLong = FRvel.Long(i);
iLat  = FRvel.Lat(i);
    arc = greatcircleArc(iLat, iLong, lat, long) * 111 ; % km
    sel = range(arc < Max_Dist);    
    sel = intersect(sel, Selected);
    if length(sel) < lim % add more stations
        add = sort(arc);
        add = add(1:lim); % add sites
        iadd = ismember(arc,add);
        inew = range(iadd);
        sel = unique(sort([sel, inew]));
        sel = intersect(sel, Selected);
    end
   V_pred_FR(i,:) = solve_LSC(iLat, iLong, lat(sel), long(sel), Vn_res(sel)*1000, Ve_res(sel)*1000,'exp1', '-v', 'bias', 'tail 0', 'no corr')'/1000;   % LSC        
end
toc

%% Interpolate LSC France velocities
clc
tic
Selected_FR = [2:14,16:75]; % Exlude CHMX, ALES, ZIMM
range = 1:length(lat);
Max_Dist = 150; % km
lim = 10;
p = 0;
step = 0.25;
for iLong = 3:step:10
    for iLat =  42:step:48
        arc = greatcircleArc(iLat, iLong, FRvel.Lat, FRvel.Long) * 111 ; % km
        sel = range(arc < Max_Dist);    
        sel = intersect(sel, Selected_FR);
        if length(sel) < lim % add more stations
            add = sort(arc);
            add = add(1:lim); % ad sites
            iadd = ismember(arc,add);
            inew = range(iadd);
            sel = unique(sort([sel, inew]));
            sel = intersect(sel, Selected_FR);
        end
        if length(sel) > 1 % do LSC !!! 
            p = p + 1; % point Number
            V_pred_grid_FR(p,:) = solve_LSC(iLat, iLong, FRvel.Lat(sel), FRvel.Long(sel), FRvel.Vn(sel), FRvel.Ve(sel),'exp1', '-v', 'bias', 'tail 0', 'no corr')'/1000;   % LSC 
            LatGrid_FR(p,1)  = iLat;
            LongGrid_FR(p,1) = iLong;
        end
    end
end
t2 = toc

%% Validation Error
FR_LSC_error = [FRvel.Vn/1000, FRvel.Ve/1000] - V_pred_FR;

close all
figure(1)
hold on; grid on
plot(FR_LSC_error(:,2), FR_LSC_error(:,1),  '.r');
plot(FRvel.Ve/1000,     FRvel.Vn/1000,      '.b');
plot(V_pred_FR(:,2),    V_pred_FR(:,1),     '.m')

%% 
% Plot Velocity field 2D map
clc
s = 500;  % [mm/yr]
try
    close (fig7)
end
fig7 = figure(7);
clr = lines(8);
hold on

% geoshow(Etopo_Europe, refvec_Etopo, 'DisplayType', 'texturemap');
% demcmap(Etopo_Europe);
% cptcmap('Europe')
Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
plot(Adriatics(:,1),Adriatics(:,2) , '--k')
% plot(FRvel.Long,      FRvel.Lat, 'or')
quiver(long(Selected),  lat(Selected),  Ve_res(Selected)*s,         Vn_res(Selected)*s,     0, 'r', 'lineWidth',1)
% quiver(long(Selected),  lat(Selected),  zeros(length(Selected),1),  Vu_res(Selected)*s,     0, 'r', 'lineWidth',1)
% quiver(FRvel.Long,      FRvel.Lat,      V_pred_FR(:,2)*s,           V_pred_FR(:,1)*s,       0, 'Color',clr(1,:), 'lineWidth',1)
% quiver(FRvel.Long,      FRvel.Lat,      zeros(76,1),                FRvel.Vu*s/1000,        0, 'k')
quiver(LongGrid_FR,     LatGrid_FR,     V_pred_grid_FR(:,2)*s,      V_pred_grid_FR(:,1)*s,  0, 'Color',clr(4,:), 'lineWidth',1)
% quiver(LongGrid,        LatGrid,        V_pred_2(:,2)*s,            V_pred_2(:,1)*s,        0, 'Color',clr(1,:), 'lineWidth',1)

% quiver(FRvel.Long,      FRvel.Lat,      FR_LSC_error(:,2)*s,        FR_LSC_error(:,1)*s,    0, 'Color',clr(3,:), 'lineWidth',1)
quiver(FRvel.Long,      FRvel.Lat,      FRvel.Ve*s/1000,            FRvel.Vn*s/1000,        0, 'k')

text(FRvel.Long,      FRvel.Lat, FRvel.Sites)
iLong = 6.75;
iLat = 45.75
ellipce_2D([km2deg(Max_Dist,6378*cosd(iLat)), Max_Dist/111], 0, [iLong, iLat], 1)
title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
xlim([2.5 10.5])
ylim([41.5 48.5])
grid on
hold off


%%

close all
figure(3)
hold on
[MapGrid_Alp, LongGrid_Alp, LatGrid_Alp] = vector2grid(V_pred_2,  LongGrid, LatGrid);
[MapGrid_Fr,  LongGrid_Fr,  LatGrid_Fr]  = vector2grid(V_pred_grid_FR, LongGrid_FR, LatGrid_FR);

imagesc(fliplr(flipud(squeeze(MapGrid_Alp(:,:,1)))))
imagesc(fliplr(flipud(squeeze(MapGrid_Fr(:,:,2)))))


