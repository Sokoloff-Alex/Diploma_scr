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
close all
fig1 = figure(1);
subplot(7,2,1) ; hold on; title('WL'); ylabel('%')
scatter(2000+AR_WL(:,1)+AR_WL(:,2)/365, AR_WL(:,4), AR_WL(:,4), AR_WL(:,3), '.'); colorbar
subplot(7,2,2); hold on; title('NL'); ylabel('%')
scatter(2000+AR_NL(:,1)+AR_NL(:,2)/365, AR_NL(:,4), AR_NL(:,4), AR_NL(:,3), '.'); colorbar

subplot(7,2,3); hold on; ; title('L5 G'); ylabel('%')
scatter(2000+AR_L5_G(:,1)+AR_L5_G(:,2)/365, AR_L5_G(:,4), AR_L5_G(:,4), AR_L5_G(:,3), '.'); colorbar
subplot(7,2,5); hold on; title('L5 R'); ylabel('%')
scatter(2000+AR_L5_R(:,1)+AR_L5_R(:,2)/365, AR_L5_R(:,4), AR_L5_R(:,4), AR_L5_R(:,3), '.'); colorbar
subplot(7,2,7); hold on; title('L5 GR'); ylabel('%')
scatter(2000+AR_L5_GR(:,1)+AR_L5_GR(:,2)/365, AR_L5_GR(:,4), AR_L5_GR(:,4), AR_L5_GR(:,3), '.'); colorbar

subplot(7,2,4); hold on; ; title('L3 G'); ylabel('%')
scatter(2000+AR_L3_G(:,1)+AR_L3_G(:,2)/365, AR_L3_G(:,4), AR_L3_G(:,4), AR_L3_G(:,3), '.'); colorbar
subplot(7,2,6); hold on; title('L3 R'); ylabel('%')
scatter(2000+AR_L3_R(:,1)+AR_L3_R(:,2)/365, AR_L3_R(:,4), AR_L3_R(:,4), AR_L3_R(:,3), '.'); colorbar
subplot(7,2,8); hold on; title('L3 GR'); ylabel('%')
scatter(2000+AR_L3_GR(:,1)+AR_L3_GR(:,2)/365, AR_L3_GR(:,4), AR_L3_GR(:,4), AR_L3_GR(:,3), '.'); colorbar

subplot(7,2,9); hold on; ; title('QIF G'); ylabel('%')
scatter(2000+AR_QIF_G(:,1)+AR_QIF_G(:,2)/365, AR_QIF_G(:,4), AR_QIF_G(:,4), AR_QIF_G(:,3), '.'); colorbar
subplot(7,2,11); hold on; title('QIF R'); ylabel('%')
scatter(2000+AR_QIF_R(:,1)+AR_QIF_R(:,2)/365, AR_QIF_R(:,4), AR_QIF_R(:,4), AR_QIF_R(:,3), '.'); colorbar
subplot(7,2,13); hold on; title('QIF GR')
scatter(2000+AR_QIF_GR(:,1)+AR_QIF_GR(:,2)/365, AR_QIF_GR(:,4), AR_QIF_GR(:,4), AR_QIF_GR(:,3), '.'); colorbar
xlabel('time, [year]')

subplot(7,2,10); hold on; ; title('L12 G'); ylabel('%')
scatter(2000+AR_L12_G(:,1)+AR_L12_G(:,2)/365, AR_L12_G(:,4), AR_L12_G(:,4), AR_L12_G(:,3), '.'); colorbar
subplot(7,2,12); hold on; title('L12 R'); ylabel('%')
scatter(2000+AR_L12_R(:,1)+AR_L12_R(:,2)/365, AR_L12_R(:,4), AR_L12_R(:,4), AR_L12_R(:,3), '.'); colorbar
subplot(7,2,14); hold on; title('L12 GR'); ylabel('%')
scatter(2000+AR_L12_GR(:,1)+AR_L12_GR(:,2)/365, AR_L12_GR(:,4), AR_L12_GR(:,4), AR_L12_GR(:,3), '.'); colorbar
xlabel('time, [year]')

print(fig1, '-depsc','-r300','AmbResStatistics.eps')
print(fig1, '-dpng','-r300','AmbResStatistics.png')






