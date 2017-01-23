% Aplen_plate_motion_est

clear all; close all; fclose all; 
clc

%%
Omega_app_Eur0     = [  0.70, -23.19, 0.038 *10^-6]';  % [deg, deg, deg/yr]' 
Omega_NNR_NUVEL_1A = [ 50.60 -112.30, 0.234 *10^-6]';  % closest !!! 
Omega_DeMet        = [ 61.06   -85.8  0.8591*10^-6]';

Re = 6371*1000;

%%
ALP_NET_CRD = readCRD('../STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('../STA/FMC_IGB_W7.VEL');

Orogen_Alp  = importOrogen('dat/PB2002_orogen_Alps.txt');

%%
CRD = cell2mat( ALP_NET_CRD(:,4:6));
VEL = cell2mat( ALP_NET_VEL(:,4:6));
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
names = ALP_NET_CRD(:,2);
flags = ALP_NET_CRD(:,7);
range_flag_W = 1:length(flags);
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);
range_flag_A = 1:length(flags);
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);
range_flag = sort([range_flag_A, range_flag_W]);

points = range_flag_W([1,4,6:end])

ALP_NET_CRD(points,:)

%% Refinment of selection,  Exclude Apls orogen
points_ref_2 = range_flag;
sites_East   = {'AXPV','PIGN','TRF2','TRFB','VARM','MEDI','GRAZ','WIEN','SPRN','HKBL','VLKM','GRS1','ZADA','CAME','UNPG','PRAT','IGMI','MOPS','GARI','PADO','VENO','VEN1','NOVE','PAZO','MDEA','JOAN','UDIN','UDI1','BOSC','MBEL','SUSE','CANV','MAVE','MPRA','CODR','PAZO','ACOM','VLCH','KOET','KOE2','ZOUF','FUSE','GRAS','VARM','AFAL','ELMO','HELM','FDOS','MPRA'};
Sites_Center = {'PF12','KTZ2','BZRG','BZ4G','BZ5G','BZ6G','RO1E','RO2E','RO3E','CO2O','ROVE','MOCA','BZRG','MITT','PAT2','PATK','HFLK','HFL2','KTZB','SBGZ','WART','HRIE','FAHR','OBE2','OBE4','VERN','KRBG','HGRA','PFAN','PFA2','BREI','PORA','SOND','COMO','CARZ','DEVE','ZIMM','ZIM2','ZI22','ZI32'};
Sites_West   = {'LI2Z','LI3Z','LINZ','WELS','RIED','ROHR','SB12','SP1N','SP2N','ZO1F','GS31','GS41','ALPE','TROP','GINA','NICA','PARO','GENO','POGG','JUVI','CHMX','FERH','ROSD','AGNE','CHTL','MODA','STEY','CHAM','JANU','PUYA','OGAG','BURE','GUIL','ACCE','CLAP','RABU','LAJA','BUIS','RSTL','MICH','BLIX','VALD','NICE','SOPH'};
Sites_Italy_and_SE         = {'TO2I','TO4I','TORI','OATO','IENG','GE1O','IE1G','VE11','TRIE','GA1I','GA2I','MO1S','MO3S','IG1I','IG2I','IG3I','UN5G','CAME','UNPG','BOLG','IGMI','PRAT','MEDI','MOPS','GARI','ZADA','PORE'};
Sites_Coast_and_Outliers   = {'GROI','MANS','BIWI','AJAC','AJ1C','AJ2C','AJ3C','CAZE','SOUR'};

Excluded = [Sites_West, sites_East, Sites_Center, Sites_Coast_and_Outliers, Sites_Italy_and_SE];

list = Excluded;
for i = 1:length(list)
    points_ref_2 = points_ref_2(strcmp(names(points_ref_2), list(i)) == 0);
end

ALP_NET_CRD(points_ref_2,:);
iiExclude = setdiff([1:length(ALP_NET_CRD)], points_ref_2);



%%
try
    close (fig1)
end
fig1 = figure(1);
scale = 1000;
hold on
grid on
Earth_coast(2)
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'-.m')
xlim([-6 18])
ylim([41 53])
plot(long(range_flag), lat(range_flag), '.b')
plot(long(points_ref_2), lat(points_ref_2), '*r')
% quiver(long(range_flag),   lat(range_flag),   Ve(range_flag)*scale,   Vn(range_flag)*scale,   0, 'b', 'lineWidth', 2)
% quiver(long(points_ref_2), lat(points_ref_2), Ve(points_ref_2)*scale, Vn(points_ref_2)*scale, 0, 'r', 'lineWidth', 2)

%% Estimate plate Euler pole

[Omega_Est, dOmega_k, DOP,Omega_Est_stack, dOmega_k_stack] = plate_motion(Omega_NNR_NUVEL_1A,'ECEF', CRD(points_ref_2,:),VEL(points_ref_2,:), 1000);

%% get residual velocities alternative

[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Est);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

%% refine selection

std_vn_res = std(Vn_res(points_ref_2),1);
std_ve_res = std(Ve_res(points_ref_2),1);
disp(['std_vn_res = ', num2str(std_vn_res*1000), ' [mm/yr]'])
disp(['std_ve_res = ', num2str(std_ve_res*1000), ' [mm/yr]'])

% in box
points_ref_2_n    = points_ref_2(abs(Vn_res(points_ref_2)) < std_vn_res); 
points_ref_2_e    = points_ref_2(abs(Ve_res(points_ref_2)) < std_ve_res);
points_ref_3 = intersect(points_ref_2_n, points_ref_2_e);

% in circle
points_ref_3= points_ref_2( sqrt( Vn_res(points_ref_2).^2 + Ve_res(points_ref_2).^2 ) < sqrt(std_vn_res^2 + std_ve_res^2))

%% Statistics
mean_en_all   = [mean(Ve_res((range_flag))), mean(Vn_res(range_flag))];
mean_en_sel_2 = [mean(Ve_res(points_ref_2)), mean(Vn_res(points_ref_2))];
mean_en_sel_3 = [mean(Ve_res(points_ref_3)), mean(Vn_res(points_ref_3))];


%%

%%
% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;
clc

% Defaults for this blog post
width =  8;     % Width in inches
height = 7.7;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw =  1;      % LineWidth
msz = 14;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

try
    close (fig9);
end
fig9 = figure(9);
subplot(4,4,[1:3,5:7,9:11])
hold on
grid on
axis equal
% axis square
pl1 = plot(Ve_res(range_flag)*1000, Vn_res(range_flag)*1000, '.b', 'MarkerSize',msz);
% text(v_long_res(range_flag)*1000, v_lat_res(range_flag)*1000, names(range_flag),'HorizontalAlignment','right')
pl2 = plot(Ve_res(points_ref_2)*1000, Vn_res(points_ref_2)*1000, '.','Color',[0 .35 0],'MarkerSize',msz);
pl3 = plot(Ve_res(points_ref_3)*1000, Vn_res(points_ref_3)*1000, '.r','MarkerSize',msz);
% plot(mean_en_all(1)*1000,   mean_en_all(2)*1000,   'xg','lineWidth',5)
plot(mean_en_sel_2(1)*1000, mean_en_sel_3(2)*1000, 'xk','MarkerSize',7)
pl4 = error_ellipse([std_ve_res^2,0;0,std_vn_res^2]*1000*1000, mean_en_sel_2*1000, 0.683,1, '--r');
pl5 = error_ellipse([std_ve_res^2,0;0,std_vn_res^2]*1000*1000, mean_en_sel_2*1000, 0.955,1, '--b');
pl6 = error_ellipse([std_ve_res^2,0;0,std_vn_res^2]*1000*1000, mean_en_sel_2*1000, 0.997,1, '--k');
ylabel('velocity SN, [mm/yr]')
title('Residual velocity')
legend([pl1 pl2 pl3 pl4 pl5 pl6], 'Alps Orogen', 'selected', 'selected: 1 sigma', '1 sigma', '2 sigma', '3 sigma')

xlim([-1 1.2])
ylim([-1 1.2])
hold off
subplot(4,4,[13:15])
hold on
[counts,bins] = hist(Ve_res(points_ref_2)*1000,[-1:0.05:1]);
bar(bins+0.02,counts,0.75,'FaceColor',[0 .5 0])
[counts,bins] = hist(Ve_res(points_ref_3)*1000,[-1:0.05:1]);
bar(bins,counts,0.75,'r')
xlabel('velocity WE, [mm/yr]')
ylabel('# of obs')
xlim([-1 1.2])
grid on
subplot(4,4,[4,8,12])
hold on
[counts,bins] = hist(Vn_res(points_ref_2)*1000,[-1:0.05:1]); %# get counts and bin locations
barh(bins+0.02,counts,0.75,'FaceColor',[0 .5 0])
[counts,bins] = hist(Vn_res(points_ref_3)*1000,[-1:0.05:1]); %# get counts and bin locations
barh(bins,counts,0.75,'r')
ylim([-1 1.2])
grid on
xlabel('# of obs')

% print(fig9, 'Residual_velocity_elipses_200.eps','-depsc','-r200');

%% Refinment of Euler pole estimation

[Omega_Est, dOmega_k2, DOP2, Omega_Est_stack2, dOmega_k_stack2] = plate_motion(Omega_Est,'ECEF',  CRD(points_ref_3,:), VEL(points_ref_3,:), 1000);
%% get residual velocities alternative

[V_res_xyz] = remove_plate_motion(CRD, VEL, Omega_Est);
[Ve_res, Vn_res, Vu_res] = XYZ2ENU(CRD,V_res_xyz); % NEU components, [m/yr m/yr m/yr]

%% plot omega
try
    close (fig3)
end
fig3 = figure(3);
subplot(2,3,[1,2,4,5])
hold on
grid on
plot(Omega_Est_stack(:,2),  Omega_Est_stack(:,1),  '.r')
plot(Omega_Est_stack2(:,2), Omega_Est_stack2(:,1), '.b')
plot(Omega_Est(2),Omega_Est(1),'om')
plot(Omega_NNR_NUVEL_1A(2),Omega_NNR_NUVEL_1A(1),'*b')
xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
title('Euler pole')
subplot(2,3,[3,6])
hold on
grid on
plot(0,Omega_NNR_NUVEL_1A(3),'*b')
plot(1:length(Omega_Est_stack), Omega_Est_stack(:,3),'.r')
plot(length(Omega_Est_stack)+1:length(Omega_Est_stack)+length(Omega_Est_stack2),Omega_Est_stack2(:,3),'.b')
title('omega')



%% MAP of Residual velocity, 2D
s = 1000; % scale
try     
    close (fig7); 
end
fig7 = figure(7);
hold on 
grid on
Earth_coast(2)
xlim([-6 18])
ylim([41 53])
plot(long(range_flag), lat(range_flag), '.b')
text(long(points_ref_3), lat(points_ref_3), names(points_ref_3))
% text(long(points_ref_3), lat(points_ref_3), names(points_ref_3),'fontsize',16)
quiver(long(range_flag),   lat(range_flag),   Ve_res(range_flag)*s,   Vn_res(range_flag)*s,   0)
quiver(long(points_ref_3), lat(points_ref_3), Ve_res(points_ref_3)*s, Vn_res(points_ref_3)*s, 0, 'r')
legend('Earth Coast','all stations','selected stations', 'Velocity','plate velocity','Residual velocity')
xlabel('Velocity EW, [mm/yr]')
ylabel('Velocity SN, [mm/yr]')
hold off

%% 3D, plot Euler pole
[xn,yn,zn] = geodetic2ecef(deg2rad(Omega_NNR_NUVEL_1A(1)),deg2rad(Omega_NNR_NUVEL_1A(2)),0,referenceEllipsoid('unitsphere'));
[xe,ye,ze] = geodetic2ecef(deg2rad(Omega_Est(1)),         deg2rad(Omega_Est(2)),         0,referenceEllipsoid('unitsphere'));


try
    close (fig2)
end
fig2 = figure(2);
scale = 100*1000*1000;
hold on
grid on
axis vis3d 
view([121 45])
Earth_coast(3)
% quiver3(CRD(range_flag,1),CRD(range_flag,2),CRD(range_flag,3), V_res_xyz(range_flag,1)*scale,V_res_xyz(range_flag,2)*scale,V_res_xyz(range_flag,3)*scale, 0)
% quiver3(CRD(points_ref_3,1),CRD(points_ref_3,2),CRD(points_ref_3,3), V_res_xyz(points_ref_3,1)*scale,V_res_xyz(points_ref_3,2)*scale,V_res_xyz(points_ref_3,3)*scale, 0)
pl1 = quiver3(CRD(points_ref_3,1),CRD(points_ref_3,2),CRD(points_ref_3,3), VEL(points_ref_3,1)*scale,VEL(points_ref_3,2)*scale,VEL(points_ref_3,3)*scale, 0, 'Linewidth',0.5);
a = 1:360;

plot3(Re*sind(a),Re*cosd(a),[zeros(length(a))],'--','Color',[.5 .5 .5])
plot3(Re*sind(a),[zeros(length(a))], Re*cosd(a),'--','Color',[.5 .5 .5])
plot3([zeros(length(a))],Re*sind(a),Re*cosd(a),'--','Color',[.5 .5 .5])

pl2 = plot3([-xn xn]*Re*1.2,[-yn yn]*Re*1.2, [-zn zn]*Re*1.2, 'k', 'lineWidth', 3 );
pl3 = plot3([-xe xe]*Re*1.2,[-ye ye]*Re*1.2, [-ze ze]*Re*1.2, 'r', 'lineWidth', 3 );


plotCircle(Omega_Est(2),Omega_Est(1),'--r')
plotCircle(Omega_NNR_NUVEL_1A(2),Omega_NNR_NUVEL_1A(1),'--k')
legend([pl1 pl2 pl3],'Station Velocities','NUVEL-1A Euler pole','Estimated Euler pole')
xlabel('x, [m]')
ylabel('y, [m]')
zlabel('z, [m]')

% print(fig2,'-dpng', '-r300','EulerPole_EU.png')



%%
for i=1:length(points_ref_3)
   fprintf('%f %f %s \n', long(points_ref_3(i)), lat(points_ref_3(i)), names{points_ref_3(i)}) 
end



