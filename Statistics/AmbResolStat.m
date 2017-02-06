% Ambiguity resolution statistics
%
% 
% load tables

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

close all
fig1 = figure(1);
subplot(4,2,1)
hold on
scatter(2000+AR_WL(:,1)+AR_WL(:,2)/365, AR_WL(:,4), AR_WL(:,4), AR_WL(:,3), 'o', 'filled')
% scatter(AR_NL(:,2), AR_NL(:,4), AR_NL(:,4), AR_NL(:,3), 'd', 'filled')
subplot(4,2,3)
plot(2000+AR_WL(:,1)+AR_WL(:,2)/365, AR_WL(:,3), 'o')

subplot(4,2,2)
hold on
scatter(AR_L5_G(:,2), AR_L5_G(:,4), AR_L5_G(:,4), AR_L5_G(:,3), 'o', 'filled')
scatter(AR_L5_R(:,2), AR_L5_R(:,4), AR_L5_R(:,4), AR_L5_R(:,3), 'd', 'filled')
% scatter(AR_L5_GR(:,2), AR_L5_GR(:,4), AR_L3_GR(:,4), AR_L3_GR(:,3), 'd', 'filled')
colorbar



