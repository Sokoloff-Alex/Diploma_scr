% Collocation
% 
% make collocation procedure to make deformation field from velocity field
% Alexandr Sokolov, KEG
% 12.10.2016

%%
clear all
close all
clc

%% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Orogen_Alp = importOrogen('dat/PB2002_orogen_Alps.txt');
lwmask = struct2array(load('dat/lwmask25.mat'));

% [Etopo_Europe, refvec_Etopo] = etopo('../../../../MAP/etopo1_bed_c_f4/etopo1_bed_c_f4.flt', 1, [40 54], [-7 19]); % ETOPO
load('ETOPO_Alps.mat')
Adriatics = struct2array(load('dat/Adriatics.mat'));
%
ALP_NET_CRD = readCRD('STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('STA/FMC_IGB_W7.VEL');

%%
flags = ALP_NET_CRD(:,7);
range_flag_W = 1:length(flags);
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);
range_flag_A = 1:length(flags);
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);
range_flag = sort([range_flag_A, range_flag_W]);

CRD = cell2mat( ALP_NET_CRD(range_flag,4:6));
VEL = cell2mat( ALP_NET_VEL(range_flag,4:6));

names = ALP_NET_CRD(range_flag,2);
DOMES = ALP_NET_CRD(range_flag,3);

% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD,VEL,names);

% remove Eurasia plate motion
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Eur);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

%% Block selection

Germany    = {'POTS','DRES','GOET','ERLA','FFMJ','WTZR','WT21','KARL', ...
              'DILL','HUEG','HOFJ'}; % Stable only
Franse_N   = {'VIGY','ERCK','BIWI','BI2I','STJ9','ENTZ','EOST', ...
              'WLBH','MOUS','AUBU','MAKS','RIXH','FRAC','LUCE','MOLV', ...
              'BUAN','ESNO','VAUC','BSCN','MLVL','MANS','HEAU','BRST', ...
              'PLOE','GROI','VFCH'}; % Stable only
Franse_S   = {'LROC','CHIZ','CAZE','TLSE','FENO','SOUR','FJCP','ARMI', ...
              'MOLA','PARD','AGDE','CAPA','PUEC','SETE','MTPL','LUVI', ...
              'BAUB','LACA','AIGL','CURA','AUMO','BANN','SAUV','STGR', ...
              'CHRN','AXPV','MARS','EGLT','SLVT','PUYV','ALLE', ...
              'TENC','ANNO','ESAB','CLFD','MONT','SAPI', ...
              'CHAR','MDOR','SJDV','MOGN','AVEN','TROC','AUTN','SIMA', ...
              'POLI','LFAZ','LEBE','STMR','PALI'}; % Stable only
Alps_SW    = {'TROP','PIGN','GINA','MICH','RSTL','BUIS','LAJA','BLIX', ...
              'VALD','SOPH','NICA','NICE','GRAS','CLAP','RABU','PQRL'};
Alps_C     = {'JUVI','JOUX','ZIMM','ZIM2','DEVE','CARZ', ...
              'COMO','PORA','BART','KRBG'};
Alps_W     = {'BURE','GUIL','JANU','PUYA','ALPE','CHAM','MODA', ...
              'CHTL','STEY','AGNE','ROSD','FCLZ','TRES'};
Alps_NE    = {'TRFB','TRF2','SPRN','RIED','LINZ', ...
              'WELS','GMND','SBG2','SBGZ','WART','HRIE','FAHR', ...
              'OBE4','BREI','HGRA','PFAN','PFA2','KTZB','KTZ2'};
Alps_E     = {'BOSC','ROVE','MOCA','FDOS','BZRG','MITT','VERN','KRGB', ...
              'HFL2','HFLK','PAT2','PATK','AFAL','ELMO','VARM'};
Austria_S  = {'KOET','KOE2','ACOM','VLCH','VLKM','GRAZ','HKBL'};
Italy_N    = {'PALO','VENE','VEN1','NOVE','PAZO','TRIE','GSR1', ...
              'MBEL','SUSE','CODR','MDEA','UDIN','UDI1','JOAN', ...
              'MPRA','FUSE','MAVE','CANV','PADO','ZOUF'};  
Apennines  = {'UNPG','IGMI','PRAT'};
Apennines_W= {'POGG','PARO','GENO','AJAC'};
Adriatic   = {'CAME','MEDI','BOLG','MOPS','ZADA','PORE','GARI'};
Adriatic_W = {'TORI','IENG','OATO'};
Outliers   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', 'ROHR'};

Stable = [Germany, Franse_N, Franse_S];

iStable      = selectRange(names, Stable);
iAlps_SW     = selectRange(names, Alps_SW);
iAlps_C      = selectRange(names, Alps_C);
iAlps_W      = selectRange(names, Alps_W);
iAlps_NE     = selectRange(names, Alps_NE);
iAlps_E      = selectRange(names, Alps_E);
iAustria_S   = selectRange(names, Austria_S);
iItaly_N     = selectRange(names, Italy_N);
iApennines   = selectRange(names, Apennines);
iApennines_W = selectRange(names, Apennines_W);
iAdriatic    = selectRange(names, Adriatic);
iAdriatic_W  = selectRange(names, Adriatic_W);

Sets = {iStable, iAlps_SW, iAlps_C, iAlps_W, iAlps_NE, iAlps_E, iAustria_S, ...
    iItaly_N, iApennines, iApennines_W, iAdriatic, iAdriatic_W};
Sets_names = {'Stable', 'Alps SW', 'Alps C', 'Alps W', 'Alps NE',  ...
    'Alps E', 'Austria S', 'Italy N', 'Apennines', 'Apennines W', 'Adriatic', 'Adriatic W'};
Selected = sort([Sets{:}]);
len = length(Sets);

%% CheckSum
% clc
disp(['# of stations selected : ',num2str(length(Selected))])
disp(['# of outliers          : ',num2str(length(Outliers))])
All = [Stable, Alps_C, Alps_E, Alps_NE, Alps_SW, Alps_W, Apennines, ...
    Apennines_W, Adriatic, Adriatic_W, Austria_S, Italy_N, Outliers];
All = sort(All)';
[iSelected, iMissing] = selectRange(names,All);
if (length(names) ~= length(iSelected)) || ~isempty(iMissing);
    length(names)
    length(iChecked)
    names(r(ismember(names,All) == 0));
    disp('Error')
else 
    disp('No stations lost')
end

%% Estimate and remove block motion from ITRF
% [Vel_res_Blocks, Vel_blocks, Omega_Blocks]= remove_block_motion(CRD, VEL, Sets);
% 
% %% Remome Plate motion from Block motion
% [V_bl_no_plate] = remove_plate_motion(CRD, Vel_blocks, Omega_Eur);
% [Ve_res_bl, Vn_res_bl, Vu_res_bl] = XYZ2ENU(CRD,V_bl_no_plate); 
% dVe = Ve_res - Ve_res_bl;
% dVn = Vn_res - Vn_res_bl;
% t = toc;
% disp(['proc time : ', num2str(t), ' [sec]'])

%% LSC, % allocation domain for grid points
clear LatGrid LongGrid V_pred V_pred*
clc
tic
range = 1:length(lat);
Max_Dist = 150; % km
lim = 10;
p = 0;
step = 0.25;
for iLong = -4:step:18
    for iLat =  42:step:53
%         if ~lwmask((90-iLat)*4+1, wrapTo360(iLong)*4+1)
            arc = greatcircleArc(iLat, iLong, lat, long) * 111 ; % km
            sel = range(arc < Max_Dist);    
            sel = intersect(sel, Selected);
            if length(sel) < lim % add more stations
                add = sort(arc);
                add = add(1:lim); % ad 2 sites, ast is prop already included
                iadd = ismember(arc,add);
                inew = range(iadd);
                sel = unique(sort([sel, inew]));
                sel = intersect(sel, Selected);
            end
            %%%%%% check to which area belong grid point %%%%%%%%%%%%%%%%%%
%             in_Adriatics = inpolygon(long(sel), lat(sel), Adriatics(:,1), Adriatics(:,2));
%             if inpolygon(iLong, iLat, Adriatics(:,1),Adriatics(:,2));
%                 sel = sel(in_Adriatics == 1);
%             else
%                 sel = sel(in_Adriatics == 0);
%             end  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if length(sel) > 1 % do LSC !!! 
                p = p + 1; % point Number
%               V_pred(p,:)   = solve_LSC_3(iLat, iLong, lat(sel), long(sel), Vn_res(sel), Ve_res(sel),  50,'no plot');   % LSC 
%               V_pred(p,:)   = solve_LSC_alt(iLat, iLong, lat(sel), long(sel), Vn_res(sel),      Ve_res(sel) ,      'no plot');   % LSC 
                V_pred_2(p,:) = solve_LSC(iLat, iLong, lat(sel), long(sel), Vn_res(sel)*1000, Ve_res(sel)*1000,'exp1', '-v', 'bias', 'tail 0', 'no corr')'/1000;   % LSC 
%                 V_pred_3(p,:) = solve_LSC(iLat, iLong, lat(sel), long(sel), Vn_res(sel)*1000, Ve_res(sel)*1000,'Hirvonen', '-v', 'no corr', 'refine')'/1000;   % LSC 
                LatGrid(p,1)  = iLat;
                LongGrid(p,1) = iLong;
            end
%         else
%             disp(['point # ', num2str(p, '%-3d'), ' :: iLat = ', num2str(iLat,'%5.2f'),' iLong = ', num2str(iLong,'%5.2f'), ' skipped, water'])
%         end
    end
end
t2 = toc


%% Plot Results 2D map
% clc
s = 250;  % [mm/yr]
try
    close (fig7)
end
clr = lines(8);
fig7 = figure(7);
hold on
xlim([-6 19])
ylim([41 53])
% xlim([2 16])
% ylim([43 49])
geoshow(Etopo_Europe, refvec_Etopo, 'DisplayType', 'texturemap');
demcmap(Etopo_Europe);
cptcmap('Europe')
% Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'--m')
plot(Adriatics(:,1),Adriatics(:,2) , '--k')
quiver(long(Selected),lat(Selected),Ve_res(Selected)*s,Vn_res(Selected)*s, 0, 'r', 'lineWidth',1)
quiver(LongGrid,      LatGrid,      V_pred_2(:,2)*s,   V_pred_2(:,1)*s,    0, 'Color',clr(1,:), 'lineWidth',1)
% ellipce_2D([km2deg(Max_Dist,6378*cosd(iLat)), Max_Dist/111], 0, [iLong, iLat], 1)
title('velocity field / deformation model') 
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Earth Coast','Alps Orogen boundary','Ardiatics','Residual velocity','LSC velocity','location','NorthWest')
% impoly(gca, [Adriatics])
grid on
hold off

%% save results
Alps_deformation = [LongGrid, LatGrid, V_pred_2(:,2), V_pred_2(:,1)];
name = 'Alps_deformation_1x1_grid';
save([name, '.mat'],'Alps_deformation');

fileID = fopen([name,'.txt'], 'w');
headString = '  Long [deg],   Lat [deg],  Vel E [m/yr],  Vel N [m/yr] \n';
formatStr = '%12.7f  %12.7f  %12.5f  %12.5f \n';
fprintf(fileID, 'Deformation / Velocity field horizontal \n');
fprintf(fileID, headString);
for i = 1:length(LongGrid)
   fprintf(fileID, formatStr, Alps_deformation(i,:)); 
end
fclose(fileID);

%% save Adriatics boundary for GMT
% %  save('Adriatics.mat' Adriatics)
% fileID = fopen('../MAP/Adriatics_boundary.txt', 'w');
% fprintf(fileID, '# Adriatics boundary MANUAL \n');
% fprintf(fileID, '# Long [deg],   Lat [deg] \n');
% for i = 1:length(Adriatics)
%    fprintf(fileID, '%12.7f  %12.7f \n',Adriatics(i,:) ); 
% end
% fclose(fileID);