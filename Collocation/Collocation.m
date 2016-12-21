% Collocation
% 
% make collocation procedure to make deformation field from velocity field
% Alexandr Sokolov, KEG
% 12.10.2016

%%
clear all
close all
clc
tic

%% Euler vector (see Alpen_plate.m)
%            lat[deg], long[deg], [deg/yr]
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
Orogen_Alp = importOrogen('bin/PB2002_orogen_Alps.txt');

%%
ALP_NET_CRD = readCRD('Velocity_field/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('Velocity_field/FMC_IGB_W7.VEL');

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

%% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD,VEL,names);

%% remove Eurasia plate motion

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
Alps_C     = {'FERR','FERH','JUVI','JOUX','ZIMM','ZIM2','DEVE','CARZ', ...
              'COMO','PORA','SOND','BART','KRBG'};
Alps_W     = {'BURE','GUIL','OGAG','JANU','PUYA','ALPE','CHAM','MODA', ...
              'CHTL','STEY','AGNE','ROSD','FCLZ','TRES'};
Alps_NE    = {'TRFB','TRF2','SPRN','ROHR','RIED','LINZ', ...
              'WELS','GMND','SBG2','SBGZ','WART','HRIE','FAHR','OBE2', ...
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
Outliers   = {'HELM', 'WIEN'};

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

Sets = {iStable, iAlps_SW, iAlps_C, iAlps_W, iAlps_NE, iAlps_E, iAustria_S, iItaly_N, iApennines, iApennines_W, iAdriatic, iAdriatic_W};
Sets_names = {'Stable', 'Alps SW', 'Alps C', 'Alps W', 'Alps NE', 'Alps E', 'Austria S', 'Italy N', 'Apennines', 'Apennines W', 'Adriatic', 'Adriatic W'};
Selected = sort([Sets{:}]);

%% CheckSum
% clc
disp(['# of stations selected : ',num2str(length(Selected))])
disp(['# of outliers          : ',num2str(length(Outliers))])
All = [Stable, Alps_C, Alps_E, Alps_NE, Alps_SW, Alps_W, Apennines, Apennines_W, Adriatic, Adriatic_W, Austria_S, Italy_N, Outliers];
All = sort(All)';
[iSelected iMissing] = selectRange(names,All);
if (length(names) ~= length(iSelected)) || ~isempty(iMissing);
    length(names)
    length(iChecked)
    names(r(ismember(names,All) == 0));
else 
    disp('All stations selected')
end

%% Test for threshold
clc
[Vel_res_Blocks, Vel_blocks, Omega_Blocks]= remove_block_motion(CRD, VEL, {iStable});

%% Estimate and remove block motion from ITRF

[Vel_res_Blocks, Vel_blocks, Omega_Blocks]= remove_block_motion(CRD, VEL, Sets);

%% Remome Plate motion from Block motion
[V_bl_no_plate] = remove_plate_motion(CRD, Vel_blocks, Omega_Eur);

[Ve_res_bl, Vn_res_bl, Vu_res_bl] = XYZ2ENU(CRD,V_bl_no_plate); 

dVe = Ve_res - Ve_res_bl;
dVn = Vn_res - Vn_res_bl;

t = toc;
disp(['proc time : ', num2str(t), ' [sec]'])

%%
% [Ve_res_bl, Vn_res_bl, Vu_res_bl] = XYZ2ENU(CRD,Vel_res_Blocks); 
% [Ve_bl,     Vn_bl,     Vu_bl]     = XYZ2ENU(CRD,Vel_blocks); 

%% MAP of Residual velocity, 2D Selection map
s = 1000; % scale
clc
len = length(Sets);
clr = cptcmap('sst','ncol', len);
try     
    close (fig1); 
end
fig1 = figure(1);
hold on 
grid on
Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'-r')
xlim([-6 18])
ylim([41 53])
title('velocity of blocks - vel of plate')
% plot(long, lat,'*r')
for i = 1:len
%     plot(  long(Sets{i}),  lat(Sets{i}),  '.','MarkerSize',8,  'Color',clr(i,:))
%     quiver(long(Sets{i}),  lat(Sets{i}),  Ve_res(Sets{i})*s,    Vn_res(Sets{i})*s,     0,'Color',clr(i,:))
    quiver(long(Sets{i}),  lat(Sets{i}),  Ve_res_bl(Sets{i})*s, Vn_res_bl(Sets{i})*s,  0,'Color',clr(i,:))
%     quiver(long(Sets{i}),  lat(Sets{i}),  dVe(Sets{i})*s,       dVn(Sets{i})*s,        0,'Color',clr(i,:))
end

% quiver(long,                lat,                Ve_res_bl*s,                Vn_res_bl*s,   0,'r')
% quiver(long,                lat,                dVe*s,                dVn*s,                0,'r')
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Earth Coast','Alps Orogen border','Residual velocity','Apls orogen','location','NorthWest')
hold off

%% Statistics for each block

dVe_mm = dVe * 1000;
dVn_mm = dVn * 1000;
for i = 1:length(Sets)
   std_dVe(i,1)  = std(dVe_mm(Sets{i}));  % std  [mm/yr] 
   std_dVn(i,1)  = std(dVn_mm(Sets{i}));  % std  [mm/yr]
   mean_dVe(i,1) = mean(dVe_mm(Sets{i})); % mean [mm/yr]
   mean_dVn(i,1) = mean(dVn_mm(Sets{i})); % mean [mm/yr]
   Corr_EN(i,1)  = corr(dVe(Sets{i}), dVn(Sets{i}));
   sigma_EN(i,1) = Corr_EN(i)*std_dVe(i)*std_dVn(i);
end

x_max = max(abs(dVe_mm(Selected)))
y_max = max(abs(dVn_mm(Selected)))
x_lim = 1;
y_lim = 1;

if x_max > 1
    x_lim = x_max + 0.2;
end
if y_max > 1
    y_lim = y_max + 0.2;
end
% plot statistics
try
    close (fig2)
end
fig2 = figure(2);
for i= 1:length(Sets)
   cov_mat = [std_dVe(i)^2 ,sigma_EN(i) ; sigma_EN(i), std_dVn(i)^2];
   mean_value = [mean_dVe(i), mean_dVn(i)];
   subplot(3,4,i)
   hold on
   grid on
   plot( dVe_mm(Sets{i}), dVn_mm(Sets{i}), '.')
   error_ellipse(cov_mat, mean_value, 0.6827, 'r') % 1 sigma, 68.27 % confidence
   error_ellipse(cov_mat, mean_value, 0.9545, 'r') % 2 sigma, 95.45 % confidence
   error_ellipse(cov_mat, mean_value, 0.9972, 'r') % 3 sigma, 99.73 % confidence
   text( dVe_mm(Sets{i}), dVn_mm(Sets{i}), names(Sets{i}))
   title([Sets_names{i}])
   xlabel('V_W_E [mm/yr]')
   ylabel('V_S_N [mm/yr]')
%    xlim([-x_lim x_lim ])
%    ylim([-y_lim y_lim ])
   hold off
end


%% build  4D feature space: lat. lon, azim, magn

[THETA,RHO] = cart2pol(Ve_res, Vn_res);

RHO_norm = RHO/max(RHO);
cm = colormap; % returns the current color map
for i = 1:length(RHO)
    colorID(i,1) = max((RHO_norm(i)),  sum( RHO_norm(i) > [0:1/length(cm(:,1)):1 ]) ); 
    myColor(i,:) = cm(colorID(i), :); % returns your color
end

try
    close (fig4)
end
fig4 = figure(4);
hold on
grid on
xlim([-6 18])
ylim([41 53])
% scatter3(long, lat, rad2deg(THETA), RHO*1000*100, myColor, 'o')
% clr = jet(12);
for i = 1:len
    scatter3(long(Sets{i}),  lat(Sets{i}), rad2deg(THETA(Sets{i})), RHO(Sets{i})*1000*10*5, repmat(clr(i,:),length(Sets{i}),1), 'filled')
    plot3(   long(Sets{i}),  lat(Sets{i}), -200*ones(size(Sets{i})), 'x', 'Color', clr(i,:) )
    plot3(Orogen_Alp(:,1),Orogen_Alp(:,2), -200*ones(size(Orogen_Alp(:,2))),'-r')
end
xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
zlabel('Azimuth, [deg]')

%% LSC, % allocation domain for grid points
clc

tic
Max_Dist = 250; % km
range = 1:length(lat);
p = 0;
step = 0.5
for iLong = -4:step:16
    for iLat = 42:step:52      
        p = p + 1; % point Number
        point = [iLong, iLat] ;
        arc = distance(point(1,2), point(1,1),lat, long) * 111 ; % km
        sel = range(arc < Max_Dist);
        sel = intersect(sel, Selected);
        d = arc(sel);
        disp(['point # ', num2str(p, '%-3d'), ' :: iLat = ', num2str(iLat,'%5.2f'),' iLong = ', num2str(iLong,'%5.2f'), ' # obs.: ', num2str(length(sel))])
        % assemble cov matrixes
        if ~isempty(sel)
            
            [C_new, C_obs] = build_LSC_CovMat_alt(iLat, iLong, lat(sel), long(sel), Vn_res(sel), Ve_res(sel), Max_Dist);
            V_obs = [Vn_res(sel); Ve_res(sel)]; % [mm/yr]  
        end
        V_pred = C_new' * C_obs^-1 * V_obs;
        Vn_pred_stack(p, 1) = V_pred(1);
        Ve_pred_stack(p, 1) = V_pred(2);
        LatGrid(p,1)  = iLat;
        LongGrid(p,1) = iLong;
    end
end

t2 = toc

%% Plot Results 2D map
clc
s = 500; 

try
    close (fig7)
end
fig7 = figure(7);
hold on
grid on
Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'-r')
xlim([-6 18])
ylim([41 53])
title('velocity of blocks - vel of plate')

% for i = 1:len
% %     plot(  long(Sets{i}),  lat(Sets{i}),  '.','MarkerSize',8,  'Color',clr(i,:))
% %     quiver(long(Sets{i}),  lat(Sets{i}),  Ve_res(Sets{i})*s,    Vn_res(Sets{i})*s,     0,'Color',clr(i,:))
% %     quiver(long(Sets{i}),  lat(Sets{i}),  Ve_res_bl(Sets{i})*s, Vn_res_bl(Sets{i})*s,  0,'Color',clr(i,:))
% %     quiver(long(Sets{i}),  lat(Sets{i}),  dVe(Sets{i})*s,       dVn(Sets{i})*s,        0,'Color',clr(i,:))
% end
% plot(long(Selected),        lat(Selected),     '.r')
% plot(LongGrid,    LatGrid, '.b')
% quiver(long(Selected),      lat(Selected),       dVe(Selected)*s,         dVn(Selected)*s, 0,'r','LineWidth',2)
quiver(LongGrid,  LatGrid, Ve_pred_stack*s, Vn_pred_stack*s,0, 'b', 'LineWidth',1)
quiver(long,                lat,                Ve_res*s,                Vn_res*s,   0,'r')
% quiver(long,                lat,                Ve_res_bl*s,                Vn_res_bl*s,   0,'r')

xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
legend('Earth Coast','Alps Orogen border','Residual velocity','LSC velocity','location','NorthWest')
hold off

print('fig7', '-dpng', '-r300');


