% Ambiguity resolution statistics
%
% 
% load tables
clc

AR_WL = read_xyz_table('../dat/OUT/AbmRes/AR_WL');
AR_NL = read_xyz_table('../dat/OUT/AbmRes/AR_NL');

AR_L5_G  = read_xyz_table('../dat/OUT/AbmRes/AR_L5_G');
AR_L5_R  = read_xyz_table('../dat/OUT/AbmRes/AR_L5_R');
AR_L5_GR = read_xyz_table('../dat/OUT/AbmRes/AR_L5_GR');

AR_L3_G  = read_xyz_table('../dat/OUT/AbmRes/AR_L3_G');
AR_L3_R  = read_xyz_table('../dat/OUT/AbmRes/AR_L3_R');
AR_L3_GR = read_xyz_table('../dat/OUT/AbmRes/AR_L3_GR');

AR_QIF_G  = read_xyz_table('../dat/OUT/AbmRes/AR_QIF_G');
AR_QIF_R  = read_xyz_table('../dat/OUT/AbmRes/AR_QIF_R');
AR_QIF_GR = read_xyz_table('../dat/OUT/AbmRes/AR_QIF_GR');

AR_L12_G  = read_xyz_table('../dat/OUT/AbmRes/AR_L12_G');
AR_L12_R  = read_xyz_table('../dat/OUT/AbmRes/AR_L12_R');
AR_L12_GR = read_xyz_table('../dat/OUT/AbmRes/AR_L12_GR');

%% plot statistics
clc

% Defaults for this blog post
width = 7;     % Width in inches
height =7;    % Height in inches
alw = 0.5;    % AxesLineWidth
fsz = 4;      % Fontsize
lw =  .05;       % LineWidth
msz = 1;       % MarkerSize

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
subplot(7,4,[1 2]) ; hold on; title('WL'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_WL(:,1)+AR_WL(:,2)/365, AR_WL(:,4), AR_WL(:,4), AR_WL(:,3), '.'); c = colorbar; 
subplot(7,4,[3 4]); hold on; title('NL');  xlim([2004 2017]);
scatter(2000+AR_NL(:,1)+AR_NL(:,2)/365, AR_NL(:,4), AR_NL(:,4), AR_NL(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');

subplot(7,4,[5 6]); hold on; ; title('L5 G'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_L5_G(:,1)+AR_L5_G(:,2)/365, AR_L5_G(:,4), AR_L5_G(:,4), AR_L5_G(:,3), '.'); c = colorbar;
subplot(7,4,[9 10]); hold on; title('L5 R'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_L5_R(:,1)+AR_L5_R(:,2)/365, AR_L5_R(:,4), AR_L5_R(:,4), AR_L5_R(:,3), '.'); c = colorbar;
subplot(7,4,[13 14]); hold on; title('L5 GR'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_L5_GR(:,1)+AR_L5_GR(:,2)/365, AR_L5_GR(:,4), AR_L5_GR(:,4), AR_L5_GR(:,3), '.'); c = colorbar; 

subplot(7,4,[7 8]); hold on; ; title('L3 G'); xlim([2004 2017]);
scatter(2000+AR_L3_G(:,1)+AR_L3_G(:,2)/365, AR_L3_G(:,4), AR_L3_G(:,4), AR_L3_G(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');
subplot(7,4,[11 12]); hold on; title('L3 R'); xlim([2004 2017]);
scatter(2000+AR_L3_R(:,1)+AR_L3_R(:,2)/365, AR_L3_R(:,4), AR_L3_R(:,4), AR_L3_R(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');
subplot(7,4,[15 16]); hold on; title('L3 GR');  xlim([2004 2017]);
scatter(2000+AR_L3_GR(:,1)+AR_L3_GR(:,2)/365, AR_L3_GR(:,4), AR_L3_GR(:,4), AR_L3_GR(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');

subplot(7,4,[17 18]); hold on; ; title('QIF G'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_QIF_G(:,1)+AR_QIF_G(:,2)/365, AR_QIF_G(:,4), AR_QIF_G(:,4), AR_QIF_G(:,3), '.'); c = colorbar; 
subplot(7,4,[21 22]); hold on; title('QIF R'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_QIF_R(:,1)+AR_QIF_R(:,2)/365, AR_QIF_R(:,4), AR_QIF_R(:,4), AR_QIF_R(:,3), '.'); c = colorbar; 
subplot(7,4,[25 26]); hold on; title('QIF GR'); ylabel('% Resolved'); xlim([2004 2017]);
scatter(2000+AR_QIF_GR(:,1)+AR_QIF_GR(:,2)/365, AR_QIF_GR(:,4), AR_QIF_GR(:,4), AR_QIF_GR(:,3), '.'); c = colorbar; 
xlabel('time, [year]')

subplot(7,4,[19 20]); hold on; ; title('L12 G');  xlim([2004 2017]);
scatter(2000+AR_L12_G(:,1)+AR_L12_G(:,2)/365, AR_L12_G(:,4), AR_L12_G(:,4), AR_L12_G(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');
subplot(7,4,[23 24]); hold on; title('L12 R');  xlim([2004 2017]);
scatter(2000+AR_L12_R(:,1)+AR_L12_R(:,2)/365, AR_L12_R(:,4), AR_L12_R(:,4), AR_L12_R(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');
subplot(7,4,[27 28]); hold on; title('L12 GR');  xlim([2004 2017]);
scatter(2000+AR_L12_GR(:,1)+AR_L12_GR(:,2)/365, AR_L12_GR(:,4), AR_L12_GR(:,4), AR_L12_GR(:,3), '.'); c = colorbar; ylabel(c,'# of Bsl.');
xlabel('time, [year]')

%%
print(fig1, '-depsc','-r300','AmbResStatistics.eps')
print(fig1, '-dpdf','-r300','AmbResStatistics.pdf')

print(fig1, '-dpng','-r300','AmbResStatistics.png')






