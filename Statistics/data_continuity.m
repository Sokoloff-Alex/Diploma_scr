% data_continuity

close all
clear all
clc

ALP_NET_CRD = readCRD('STA/FMC_IGB_W7.CRD');
ALP_NET_VEL = readVEL('STA/FMC_IGB_W7.VEL');


%
flags = ALP_NET_CRD(:,7);
range_flag_W = 1:length(flags);
range_flag_W = range_flag_W(strcmp(flags, 'W') == 1);
range_flag_A = 1:length(flags);
range_flag_A = range_flag_A(strcmp(flags, 'A') == 1);
range_flag = sort([range_flag_A, range_flag_W]);

CRD_all = cell2mat( ALP_NET_CRD(range_flag,4:6));
VEL_all = cell2mat( ALP_NET_VEL(range_flag,4:6));

names_all = ALP_NET_CRD(range_flag,2);
DOMES = ALP_NET_CRD(range_flag,3);

[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);

% Megre artificial stations
[CRD,VEL,names] = merge_stations(CRD_all,VEL_all,names_all);

% remove Eurasia plate motion
[Ve,Vn, Vu, lat, long,  h]  = XYZ2ENU(CRD,VEL);
[Ve_all,Vn_all, Vu_all, lat_all, long_all,  h_all]  = XYZ2ENU(CRD_all,VEL_all);


load('../dat/SNX/SINEX.mat')
SNX_cov = SINEX.SOLUTION.COVA_ESTIMATE;
[CovVenuSNX, SigmaVenu, CorrVen, AngleV] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'VEL');
[CovRenuSNX, SigmaRenu, CorrRen, AngleR] = SNX_cov_transformXYZ2ENU(SNX_cov,lat_all, long_all, 'CRD');
[CovVenu] = megreCov(CovVenuSNX, names_all);

[CRD, SigmaVenu_merged_old, name ] = merge_stations(CRD_all,SigmaVenu,names_all);
[CRD, SigmaVenu_merged, name ] = merge_stations_precision(CRD_all,SigmaVenu,names_all);
[CRD, AngleV_merged] = merge_stations(CRD_all,AngleV,names_all);

[CRD, SigmaRenu_merged, name ] = merge_stations_precision(CRD_all,SigmaRenu,names_all);
[CRD, AngleR_merged] = merge_stations(CRD_all,AngleR,names_all);
Angle_v = AngleV_merged(:,1);



%% compute common observation period
% 
t_start = SINEX.SOLUTION.EPOCHS.DATA_START;
t_end   = SINEX.SOLUTION.EPOCHS.DATA_END;

Sig_Vel_xyz = SINEX.SOLUTION.ESTIMATE.Data.VEL_STD;

t_start = [str2num(t_start(:,1:2)) + str2num(t_start(:,4:6))/365.25]
t_end   = [str2num(t_end(  :,1:2)) + str2num(t_end(  :,4:6))/365.25]

dt = t_end - t_start;
% 
% %%
[dtobs_sum] = merge_stations_sum(dt,names_all);

%% scatter of presicion vs tobs

% Defaults for this blog post
width =  5.5;  % Width in inches
height = 4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 7;       % Fontsize
lw =  1;       % LineWidth
msz = 14;      % MarkerSize

% The properties we've been using in the figures
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

close all

fig2 = figure(2);
subplot(2,1,1)
hold on
grid on
title('Formal Error of velocities for artificial stations')
plot(dt, SigmaVenu(:,1)*1000*sqrt(1.8)*20,'.b')
plot(dt, SigmaVenu(:,2)*1000*sqrt(1.8)*20,'.g')
plot(dt, SigmaVenu(:,3)*1000*sqrt(1.8)*20,'.r')
l1 = legend('\sigma_V_e','\sigma_V_n','\sigma_V_u')
% set(l1, 'MarkerSize',11)
ylim([0 1.6])
xlim([0 14])
set(gca,'Ytick',0:0.2:2)
set(gca,'YTickLabel',{num2str([0:0.2:2]')})
set(gca,'Xtick',0:13)
set(gca,'XTickLabel',{num2str([0:13]')})

% xlabel('observation time, [year]')
ylabel(' precision, [mm/yr]')

subplot(2,1,2)
hold on
grid on
title('Averaged formal error of station velocities for original stations')
plot(dtobs_sum, SigmaVenu_merged(:,1)*1000*sqrt(1.8)*20,'.b')
plot(dtobs_sum, SigmaVenu_merged(:,2)*1000*sqrt(1.8)*20,'.g')
plot(dtobs_sum, SigmaVenu_merged(:,3)*1000*sqrt(1.8)*20,'.r')
% legend('\sigma_V_e','\sigma_V_n','\sigma_V_u')
ylim([0 1.6])
xlim([0 14])
set(gca,'Ytick',0:0.2:2)
set(gca,'YTickLabel',{num2str([0:0.2:2]')})
set(gca,'Xtick',0:13)
set(gca,'XTickLabel',{num2str([0:13]')})
xlabel('observation time, [year]')
ylabel(' precision, [mm/yr]')

%%
print(fig2, '-depsc','-r300','VelPrecisionDistribution.eps')
print(fig2, '-dpdf','-r300','VelPrecisionDistribution.pdf')


%% plot histogramm

[nelements1,centers1] = hist(dtobs_sum, [2:1:12]);
[nelements2,centers2] = hist(dt,        [0:1:12]);

% Defaults for this blog post
width =  5.5;  % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 7;      % Fontsize
lw =  1;       % LineWidth
msz = 14;      % MarkerSize

% The properties we've been using in the figures
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


close all
fig1 = figure(1);
hold on
b1 = bar(centers1+0.2, nelements1, 0.4, 'b');
b2 = bar(centers2-0.2, nelements2, 0.4, 'g');
text(centers1, nelements1+2, num2str(nelements1'),'Color','b')
text(centers2-0.55, nelements2+2,  num2str(nelements2'),'Color',[0 .6 0])
legend([b1 b2], {['# of original stations, (total: ',num2str(length(dtobs_sum)),')'], ...
    ['# of artificial stations, (total: ',num2str(length(dt)),')']},'location','NorthEast')
xlabel('Duration of observations, years')
ylabel('# of stations')
xlim([-1 13])
ylim([0 55])
set(gca,'Xtick',0:12)
set(gca,'XTickLabel',{'>0', num2str([1:12]')})

%%
print(fig1,'-depsc','-r300', 'HistContinuity.eps')
% print(fig1,'-dpdf', '-r300', 'HistContinuity.pdf')
% print(fig1,'-dpng', '-r300', 'HistContinuity.png')


