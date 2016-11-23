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

%% Resove large vel discripancy at ZIMM ZIM2

Selected_FR = [2:14,16:75]; % Exlude CHMX, ALES, ZIMM
outFR = [1, 15, 75, 76];
FRvel.Ve(77) = ( FRvel.Ve(75) + FRvel.Ve(76) )/2;
FRvel.Vn(77) = ( FRvel.Vn(75) + FRvel.Vn(76) )/2;
FRvel.Long(77) = ( FRvel.Long(75) + FRvel.Long(76) )/2;
FRvel.Lat(77)  = ( FRvel.Lat(75)  + FRvel.Lat(76)  )/2;
FRvel.Sites(77) = {'ZIM'};


%% Interpolate LSC France velocities
clc
tic
range = 1:length(lat);
Max_Dist = 150; % km
lim = 10;
p = 0;
step = 0.25;
for iLong = 3:step:11
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
sc1 = 1000*0.5;  % [mm/yr] > [m/yr]
scAlp = 10;
scFR = sc1*100;
scSigVe = sc1 * scAlp * 100;
scSigVn = sc1 * scAlp * 350;

try
    close (fig7)
end
fig7 = figure(7);
clr = lines(8);
hold on

geoshow(Etopo_Europe, refvec_Etopo, 'DisplayType', 'texturemap');
demcmap(Etopo_Europe);
cptcmap('Europe')
% Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
plot(Adriatics(:,1),Adriatics(:,2) , '--k')
% plot(FRvel.Long,      FRvel.Lat, 'or')
quiver(wrapTo180(long(Selected)),  lat(Selected),  Ve_res(Selected)*sc1,         Vn_res(Selected)*sc1,     0, 'b', 'lineWidth',2)
% quiver(long(Selected),  lat(Selected),  zeros(length(Selected),1),  Vu_res(Selected)*s,     0, 'r', 'lineWidth',1)
quiver(LongGrid,        LatGrid,        V_pred_2(:,2)*sc1,            V_pred_2(:,1)*sc1,        0, 'b', 'lineWidth',1)

for i = 1:length(Selected) 
   ellipce_2D([km2deg(SigmaVe(Selected(i))*scSigVe, 6378*cosd(lat(Selected(i)))), SigmaVn(Selected(i))*scSigVn/111], ...
       0, [long(Selected(i)) + Ve_res(Selected(i))*sc1, lat(Selected(i))  + Vn_res(Selected(i))*sc1], 1, 'b') 
end

% quiver(FRvel.Long,      FRvel.Lat,      V_pred_FR(:,2)*s,           V_pred_FR(:,1)*s,       0, 'm', 'lineWidth',2)
% quiver(FRvel.Long,      FRvel.Lat,      zeros(76,1),                FRvel.Vu*s/1000,        0, 'k')
quiver(LongGrid_FR,     LatGrid_FR,     V_pred_grid_FR(:,2)*sc1,      V_pred_grid_FR(:,1)*sc1,  0, 'r', 'lineWidth',1)

% quiver(FRvel.Long,      FRvel.Lat,      FR_LSC_error(:,2)*s,        FR_LSC_error(:,1)*s,    0, 'Color',clr(3,:), 'lineWidth',1)
quiver(FRvel.Long,      FRvel.Lat,      FRvel.Ve*sc1/1000,            FRvel.Vn*sc1/1000,        0, 'r', 'lineWidth',2)
% quiver(FRvel.Long(outFR),FRvel.Lat((outFR)),FRvel.Ve((outFR))*s/1000,FRvel.Vn((outFR))*s/1000,        0, 'Color',[.5 .5 .5], 'lineWidth',2)

for i = 1:length(FRvel.Lat)
   ellipce_2D([ km2deg(FRvel.SigmaVe(i), 6378*cosd(FRvel.Lat(i))), km2deg(FRvel.SigmaVn(i)) ]*scFR/1000,  ...
       0, [wrapTo180(FRvel.Long(i)) + FRvel.Ve(i)*sc1/1000, FRvel.Lat(i)  + FRvel.Vn(i)*sc1/1000], 1, 'r') 
end

% text(long, lat, names)
% text(FRvel.Long,      FRvel.Lat, FRvel.Sites)
iLong = 6.75;
iLat = 45.75
% ellipce_2D([km2deg(Max_Dist,6378*cosd(iLat)), Max_Dist/111], 0, [iLong, iLat], 1)
title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
xlim([3 11])
ylim([42 48])
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


%% save results France Deformation
clc
France_deformation = [LongGrid_FR, LatGrid_FR, V_pred_grid_FR(:,2), V_pred_grid_FR(:,1)];
name = 'France_deformation_0.25x0.25_grid';
save([name, '.mat'],'France_deformation');

% write txt file
fileID = fopen([name,'.txt'], 'w');
headString = '  Long [deg],   Lat [deg],  Vel E [m/yr],  Vel N [m/yr],   SigmaVe[m/yr],    SigmaVn[m/yr] \n';
formatStr = '%12.7f  %12.7f  %12.5f  %12.5f   %14.8f   %14.8f \n';
fprintf(fileID, 'France, Deformation / Velocity field horizontal \n');
fprintf(fileID, headString);
for i = 1:size(France_deformation,1)
   fprintf(fileID, formatStr, France_deformation(i,:)); 
end
fclose(fileID);

%% save results France Velocity field
clc
France_velocity_field = [FRvel.Long, FRvel.Lat, FRvel.Ve/1000, FRvel.Vn/1000, FRvel.SigmaVe/1000 FRvel.SigmaVn/1000];
name = 'France_velocity_field';
save([name, '.mat'],'France_velocity_field');

% write txt file
fileID = fopen([name,'.txt'], 'w');
headString = '  Long [deg],   Lat [deg],  Vel E [m/yr],  Vel N [m/yr],    SigmaVe[m/yr],    SigmaVn[m/yr],   SITE \n';
formatStr = '%12.7f  %12.7f  %12.5f  %12.5f    %14.8f  %14.8f       %4s \n';
fprintf(fileID, 'France, Velocity field horizontal \n');
fprintf(fileID, headString);
for i = 1:size(France_velocity_field,1)
   fprintf(fileID, formatStr, France_velocity_field(i,:), cell2mat( FRvel.Sites(i))); 
end
fclose(fileID);


