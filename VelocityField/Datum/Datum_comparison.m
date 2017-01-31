% Datum_comparison
%
% compare CRD and VEL between datums: 
% EPN EPN_A_IGb08.SNX
% IGS IGS16P26.SNX
% IGb IGb08.SNX
%
% for sations: GRAS, GRAZ, LROC, MEDI, WTZR, ZIMM

close all
clear all
clc
%%

IGB_SNX = readSNX('../../dat/SNX/IGb08_noCov.snx', 'no COVA');
IGS_SNX = readSNX('../../dat/SNX/IGS16P26.snx', 'no COVA');
EPN_SNX = readSNX('../../dat/SNX/EPN_A_IGb08_no_COVA.SNX', 'no COVA');

%%
% save('../../dat/SNX/IGb08_noCov.SNX.mat',  'IGB_SNX')
% save('../../dat/SNX/IGS16P26.SNX.mat',     'IGS_SNX')
% save('../../dat/SNX/IGb08_no_COVA.SNX.mat','EPN_SNX')

load('../../dat/SNX/IGb08_noCov.SNX.mat')
load('../../dat/SNX/IGS16P26.SNX.mat')
load('../../dat/SNX/IGb08_no_COVA.SNX.mat')

%% compare CRD and VEL of stations: GRAS, GRAZ, LROC, MEDI, WTZR, ZIMM
set_ref_sites = {'GRAS', 'GRAZ', 'LROC', 'MEDI', 'WTZR', 'ZIMM'};
range = [1:length(IGB_SNX.SOLUTION.EPOCHS.CODE)]';
iiIGB = range(ismember(IGB_SNX.SOLUTION.EPOCHS.CODE, set_ref_sites))
range = [1:length(IGS_SNX.SOLUTION.EPOCHS.CODE)]';
iiIGS = range(ismember(IGS_SNX.SOLUTION.EPOCHS.CODE, set_ref_sites))
range = [1:length(EPN_SNX.SOLUTION.EPOCHS.CODE)]';
iiEPN = range(ismember(EPN_SNX.SOLUTION.EPOCHS.CODE, set_ref_sites))

IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:)
IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:)
EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:)

%% sort 

iiIGS = [iiIGS([10:12, 1:5, 15, 13:14, 8:9, 6:7 ])]'
iiEPN = [iiEPN([1,8:9,4,6,10:12,7,2,13,3,5])],

[IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:) , IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:)]
IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:)
IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:)
EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:)

%%

IGB_crd_set = [IGB_SNX.SOLUTION.ESTIMATE.Data.CRD(iiIGB,:)]
IGB_vel_set = [IGB_SNX.SOLUTION.ESTIMATE.Data.VEL(iiIGB,:)]

IGS_crd_set = [IGS_SNX.SOLUTION.ESTIMATE.Data.CRD(iiIGS,:)]
IGS_vel_set = [IGS_SNX.SOLUTION.ESTIMATE.Data.VEL(iiIGS,:)]

EPN_crd_set = [EPN_SNX.SOLUTION.ESTIMATE.Data.CRD(iiEPN,:)]
EPN_vel_set = [EPN_SNX.SOLUTION.ESTIMATE.Data.VEL(iiEPN,:)]

%% average

IGB_crd_set = [mean(IGB_crd_set(1:3,:))
               mean(IGB_crd_set(4:8,:))
               IGB_crd_set(9,:)
               mean(IGB_crd_set(10:11,:))
               mean(IGB_crd_set(12:13,:))
               mean(IGB_crd_set(14:15,:))];

IGB_vel_set = [mean(IGB_vel_set(1:3,:))
               mean(IGB_vel_set(4:8,:))
               IGB_vel_set(9,:)
               mean(IGB_vel_set(10:11,:))
               mean(IGB_vel_set(12:13,:))
               mean(IGB_vel_set(14:15,:))];
           
IGS_crd_set = [mean(IGS_crd_set(1:3,:))
               mean(IGS_crd_set(4:8,:))
               IGS_crd_set(9,:)
               mean(IGS_crd_set(10:11,:))
               mean(IGS_crd_set(12:13,:))
               mean(IGS_crd_set(14:15,:))];

IGS_vel_set = [mean(IGS_vel_set(1:3,:))
               mean(IGS_vel_set(4:8,:))
               IGS_vel_set(9,:)
               mean(IGS_vel_set(10:11,:))
               mean(IGS_vel_set(12:13,:))
               mean(IGS_vel_set(14:15,:))];

EPN_crd_set = [mean(EPN_crd_set(1:3,:))
               mean(EPN_crd_set(4:8,:))
               EPN_crd_set(9,:)
               mean(EPN_crd_set(10:11,:))
               mean(EPN_crd_set(12:13,:))];

EPN_vel_set = [mean(EPN_vel_set(1:3,:))
               mean(EPN_vel_set(4:8,:))
               EPN_vel_set(9,:)
               mean(EPN_vel_set(10:11,:))
               mean(EPN_vel_set(12:13,:))];

%% convert xyz -> enu

[IGB_e, IGB_n, IGB_u, IGB_lat, IGB_long,  IGB_h]  = XYZ2ENU(IGB_crd_set,IGB_vel_set*1000);
[IGS_e, IGS_n, IGS_u, IGS_lat, IGS_long,  IGS_h]  = XYZ2ENU(IGS_crd_set,IGS_vel_set*1000);
[EPN_e, EPN_n, EPN_u, EPN_lat, EPN_long,  EPN_h]  = XYZ2ENU(EPN_crd_set,EPN_vel_set*1000);

%% plot 
close all
clc

%
% Defaults for this blog post
width = 8.5;     % Width in inches
height = 7;    % Height in inches
fsz = 10;      % Fontsize
lw =  1;       % LineWidth
msz = 4;       % MarkerSize

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

% close all

clr = lines(6);
ii = [1:3,5,6];

fig1 = figure(1);
subplot(3,3,1)
hold on; grid on; axis equal; axis square
de = -[IGB_e([1:3,5,6]) - EPN_e];
xlim([17.5 22.5])
ylim([17.5 22.5])
plot([17.5 22.5], [17.5 22.5])
for i = 1:5
    plot(IGB_e(ii(i)), EPN_e(i),'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de),'%5.3f'), ' mm/yr'])

subplot(3,3,2)
hold on; grid on; axis equal; axis square
dn = -[IGB_n([1:3,5,6]) - EPN_n];
xlim([15 17])
ylim([15 17])
plot([15 17], [15 17])
for i = 1:5
    plot(IGB_n(ii(i)), EPN_n(i),'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.3f'), ' mm/yr'])

subplot(3,3,3)
hold on; grid on; axis equal; axis square
du = -[IGB_u([1:3,5,6]) - EPN_u];
xlim([-.5 1.5])
ylim([-.5 1.5])
plot([-.5 1.5], [-.5 1.5])
for i = 1:5
    plot(IGB_u(ii(i)), EPN_u(i),'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(du),'%5.3f'), ' mm/yr'])

subplot(3,3,4)
hold on; grid on; axis equal; axis square
de = [IGB_e - IGS_e];
for i = 1:6
    plot(IGB_e(i), IGS_e(i),'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:))
end
xlabel('IGB')
ylabel('IGS')
xlim([18 23])
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de),'%5.2f'), ' mm/yr'])

subplot(3,3,5)
hold on; grid on; axis equal; axis square
dn = [IGB_n - IGS_n];
for i = 1:6
    plot(IGB_n(i), IGS_n(i),'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:))
end
xlabel('IGB')
ylabel('IGS')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.2f'), ' mm/yr'])
xlim([15.5 18])

s6 = subplot(3,3,6);
hold on; grid on; axis equal; axis square
du = [IGB_u - IGS_u];
for i = 1:6
    plt(i) = plot(IGB_u(i), IGS_u(i),'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:));
end
xlabel('IGB')
ylabel('IGS')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(du),'%5.2f'), ' mm/yr'])
xlim([-2 1.5])

subplot(3,3,7)
hold on; grid on; axis equal; axis square
de = [IGS_e([1:3,5,6]) - EPN_e];
for i = 1:5
    plot(EPN_e(i), IGS_e(ii(i)), 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de),'%5.2f'), ' mm/yr'])
xlim([18 23])

subplot(3,3,8)
hold on; grid on; axis equal; axis square
dn = [IGS_n([1:3,5,6]) - EPN_n];
for i = 1:5
    plot(EPN_n(i), IGS_n(ii(i)), 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.2f'), ' mm/yr'])
xlim([15.5 18])

s9 = subplot(3,3,9);
hold on; grid on; axis equal; axis square
du = [IGS_u([1:3,5,6]) - EPN_u];
for i = 1:5
    plot(EPN_u(i), IGS_u(ii(i)), 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(du),'%5.2f'), ' mm/yr'])
leg = legend(s6, [plt], 'GRAS','GRAZ','LROC','MEDI','WTZR','ZIMM','location','NorthWest');
set(leg,'FontSize',6)
xlim([-2 1.5])

%
print(fig1,'-depsc', '-r300','ComparisonOfRefSites_vel.eps')


%% diff

[IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB([1:9,12:15]),:), EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:)]
[IGB_SNX.SOLUTION.EPOCHS.DATA_START(iiIGB([1:9,12:15]),:), EPN_SNX.SOLUTION.EPOCHS.DATA_START(iiEPN,:)]
[IGB_SNX.SOLUTION.EPOCHS.DATA_END(iiIGB([1:9,12:15]),:), EPN_SNX.SOLUTION.EPOCHS.DATA_END(iiEPN,:)]


rms([IGB_vel_set - IGS_vel_set])*1000
rms([IGB_vel_set([1:9,12:15],:) - EPN_vel_set])*1000


close all
figure(1)
subplot(3,3,1)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set(:,1) , IGS_vel_set(:,1)]))*1000;
maxA = max(max([IGB_vel_set(:,1) , IGS_vel_set(:,1)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set(:,1)*1000,   IGS_vel_set(:,1)*1000, 'x' )
text(IGB_vel_set(:,1)*1000+3, IGS_vel_set(:,1)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:),'Color', 'b')
% text(IGB_vel_set(:,1)*1000-2, IGS_vel_set(:,1)*1000, IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:),'Color', 'r')
title(['rms_x = ', num2str(mean( IGB_vel_set(:,1) - IGS_vel_set(:,1))*1000) ])


subplot(3,3,4)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set(:,2) , IGS_vel_set(:,2)]))*1000;
maxA = max(max([IGB_vel_set(:,2) , IGS_vel_set(:,2)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set(:,2)*1000,   IGS_vel_set(:,2)*1000, 'x' )
text(IGB_vel_set(:,2)*1000+3, IGS_vel_set(:,2)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:),'Color', 'b')
% text(IGB_vel_set(:,2)*1000-2, IGS_vel_set(:,2)*1000, IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:),'Color', 'r')
title(['rms_x = ', num2str(mean( IGB_vel_set(:,2) - IGS_vel_set(:,2))*1000) ])

subplot(3,3,7)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set(:,3) , IGS_vel_set(:,3)]))*1000;
maxA = max(max([IGB_vel_set(:,3) , IGS_vel_set(:,3)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set(:,3)*1000,   IGS_vel_set(:,3)*1000, 'x' )
text(IGB_vel_set(:,3)*1000+3, IGS_vel_set(:,3)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB,:),'Color', 'b')
% text(IGB_vel_set(:,3)*1000-2, IGS_vel_set(:,3)*1000, IGS_SNX.SOLUTION.EPOCHS.CODE(iiIGS,:),'Color', 'r')
title(['rms_x = ', num2str(mean( IGB_vel_set(:,3) - IGS_vel_set(:,3))*1000) ])


subplot(3,3,2)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set([1:9,12:15],1) , EPN_vel_set(:,1)]))*1000;
maxA = max(max([IGB_vel_set([1:9,12:15],1) , EPN_vel_set(:,1)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set([1:9,12:15],1)*1000,   EPN_vel_set(:,1)*1000, 'x' )
text(IGB_vel_set([1:9,12:15],1)*1000+1, EPN_vel_set(:,1)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB([1:9,12:15]),:),'Color', 'b')
text(IGB_vel_set([1:9,12:15],1)*1000-1, EPN_vel_set(:,1)*1000, EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:),             'Color', 'r')
title(['rms_x = ', num2str(mean( IGB_vel_set([1:9,12:15],1) - EPN_vel_set(:,1))*1000) ])

subplot(3,3,5)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set([1:9,12:15],2) , EPN_vel_set(:,2)]))*1000;
maxA = max(max([IGB_vel_set([1:9,12:15],2) , EPN_vel_set(:,2)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set([1:9,12:15],2)*1000,     EPN_vel_set(:,2)*1000, 'x' )
% text(IGB_vel_set([1:9,12:15],2)*1000+0.2, EPN_vel_set(:,2)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB([1:9,12:15]),:),'Color', 'b')
% text(IGB_vel_set([1:9,12:15],2)*1000-0.2, EPN_vel_set(:,2)*1000, EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:),             'Color', 'r')
title(['rms_y = ', num2str(mean( IGB_vel_set([1:9,12:15],2) - EPN_vel_set(:,2))*1000) ])

subplot(3,3,8)
hold on
grid on
axis equal
minA = min(min([IGB_vel_set([1:9,12:15],3) , EPN_vel_set(:,3)]))*1000;
maxA = max(max([IGB_vel_set([1:9,12:15],3) , EPN_vel_set(:,3)]))*1000;
plot([minA maxA], [minA maxA])
plot(IGB_vel_set([1:9,12:15],3)*1000,     EPN_vel_set(:,3)*1000, 'x' )
% text(IGB_vel_set([1:9,12:15],3)*1000+0.2, EPN_vel_set(:,3)*1000, IGB_SNX.SOLUTION.EPOCHS.CODE(iiIGB([1:9,12:15]),:),'Color', 'b')
% text(IGB_vel_set([1:9,12:15],3)*1000-0.2, EPN_vel_set(:,3)*1000, EPN_SNX.SOLUTION.EPOCHS.CODE(iiEPN,:),             'Color', 'r')
title(['rms_z = ', num2str(mean( IGB_vel_set([1:9,12:15],3) - EPN_vel_set(:,3))*1000) ])

