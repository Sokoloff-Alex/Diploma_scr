% Oscillations

% FD_ALP_1.EVT
PERI365.Properties.VariableNames = {'Name','code','flag','event','Period','North','East','Up','S','remark'};


% save('PERI365','PERI365')
% save('PERP365','PERP365')


% load('PERI365')
% load('PERP365')

%%

names_t = cell2mat(PERI365.Name);

r = [1:length(names_t)]';
ii = r(ismember(names_t, names))
%%
close all
figure(5)
subplot(1,3,1)
hold on
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])
scatter(wrapTo180(long), lat, PERI365.East(ii)*10000+1, PERI365.East(ii)*1000 ,'filled' )
title('East')
c = colorbar;
ylabel(c,'mm')
caxis([0 5])

subplot(1,3,2)
hold on
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])
scatter(wrapTo180(long), lat, PERI365.North(ii)*10000+1, PERI365.North(ii)*1000 ,'filled' )
title('North')
% text(wrapTo180(long), lat, names)
c = colorbar;
ylabel(c,'mm')
caxis([0 5])

subplot(1,3,3)
hold on
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])
scatter(wrapTo180(long), lat, PERI365.Up(ii)*10000+1, PERI365.Up(ii)*1000 ,'filled' )
title('Up')
c = colorbar;
ylabel(c,'mm')
caxis([0 5])




%%
clc
close all

% Defaults for this blog post
width = 10;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw =  1;       % LineWidth
msz = 6;       % MarkerSize

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


% %% mag + phase
% clc
% close all
% fig2 = figure(2);
% colormap('hsv')
% subplot(1,3,1); hold on
% scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName8*10000+1, PERP365.VarName8 ,'filled' )
% title('East')
% xlim([-4.5 17])
% subplot(1,3,2); hold on
% scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName9*10000+1, PERP365.VarName9 ,'filled' )
% title('North')
% xlim([-4.5 17])
% subplot(1,3,3); hold on
% scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName10*10000+1, PERP365.VarName10 ,'filled' )
% title('Up')
% xlim([-4.5 17])

%
close all
fig3 = figure(3);
colormap('hsv')
subplot(1,10,[1 2 3]); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.East*10000+1, PERI365.East*1000 ,'filled' )
title('East')
% c = colorbar;
% ylabel(c,'mm')
caxis([0 5])
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])
subplot(1,10,[4 5 6]); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.North*10000+1, PERI365.North*1000 ,'filled' )
title('North')
% c = colorbar;
% ylabel(c,'mm')
caxis([0 5])
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])
subplot(1,10,[ 7 8 9 10 ]); hold on
%%
close all
figure(4); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.Up*10000+1, PERI365.Up*1000 ,'filled' )
text(wrapTo180(long_all(ii)), lat_all(ii)+0.2, names_all(ii))
text(wrapTo180(long_all(ii)), lat_all(ii), num2str(PERI365.Up*1000, '%3.1f'))
title('Up')
s = colorbar;
ylabel(s,'mm')
caxis([0 5])
Earth_coast(2)
xlim([-4.5 17])
ylim([41.8 52])

%%
close all
fig3 = figure(3);
colormap('hsv')
subplot(2,9,[10 11 ]); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), 5*ones(274,1), PERP365.East ,'filled' )
title('East')
Earth_coast(2)
xlim([-4.5 17])
ylim([41 53])
subplot(2,9, [12 13]); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), 5*ones(274,1), PERP365.North ,'filled' )
title('North')
xlim([-4.5 17])
subplot(2,9 ,[14 15 16 ]); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), 5*ones(274,1), PERP365.Up ,'filled' )
title('Up')
xlim([-4.5 17])
colorbar

%%

close all
fig3 = figure(3);
colormap('hsv')
subplot(1,4,1); hold on
s = 40;
scatter(wrapTo180(long_all(ii)), lat_all(ii), s*ones(274,1), PERI365.East*1000 ,'filled' )
title('East')
Earth_coast(2)
xlim([-4.5 17])
ylim([41 53])
subplot(1,4,2 ); hold on
Earth_coast(2)
xlim([-4.5 17])
ylim([41 53])
scatter(wrapTo180(long_all(ii)), lat_all(ii), s*ones(274,1), PERI365.North*1000 ,'filled' )
title('North')
xlim([-4.5 17])
subplot(1,4,[3 4]); hold on
Earth_coast(2)
xlim([-4.5 17])
ylim([41 53])
scatter(wrapTo180(long_all(ii)), lat_all(ii), s*ones(274,1), PERI365.Up*1000 ,'filled' )
title('Up')
xlim([-4.5 17])
colorbar

%%
print(fig3, 'dat/Pics/AnnualOscillation.eps','-depsc','-r300');
print(fig3, 'dat/Pics/AnnualOscillation.pdf','-dpdf','-r300');

