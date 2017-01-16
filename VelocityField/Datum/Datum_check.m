% Datum check

% save('../../dat/DatumcheckFMCIGBW7');
load('../../dat/DatumcheckFMCIGBW7');
Datum_dr = DatumcheckFMCIGBW7;

%%
clc
[XYZ] = lla2ecef([45, 8, 0]);
[de, dn, du] = XYZ2ENU(repmat(XYZ,length(Datum_dr),1), Datum_dr(:,3:5));

%%

%%
% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;
clc

% Defaults for this blog post
width = 7.5;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw =  1;      % LineWidth
msz = 6;       % MarkerSize

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

%
close all
fig1 = figure(1);
hold on
grid on
% plot(Datum_dr(:,1),Datum_dr(:,2)*1000,'.k')
plot(Datum_dr(:,1),Datum_dr(:,3)*1000+10,'.r')
plot(Datum_dr(:,1),Datum_dr(:,4)*1000+20,'.','Color',[0 .5 0])
plot(Datum_dr(:,1),Datum_dr(:,5)*1000+30,'.b')
plot(Datum_dr(:,1), de*1000+40, '.r')
plot(Datum_dr(:,1), dn*1000+50, '.','Color',[0 .5 0])
plot(Datum_dr(:,1), du*1000+60, '.b')
% text(660,2,   'RMS')
text(660,1+10, ['dX, rms = ', num2str( rms(Datum_dr(:,3))*1000,'%3.2f'), ' mm'])
text(660,1+20, ['dY, rms = ', num2str( rms(Datum_dr(:,4))*1000,'%3.2f'), ' mm'])
text(660,1+30, ['dZ, rms = ', num2str( rms(Datum_dr(:,5))*1000,'%3.2f'), ' mm'])
text(660,1+40, ['dE, rms = ', num2str( rms(de)*1000,'%3.2f'), ' mm'])
text(660,1+50, ['dN, rms = ', num2str( rms(dn)*1000,'%3.2f'), ' mm'])
text(660,1+60, ['dU, rms = ', num2str( rms(du)*1000,'%3.2f'), ' mm'])
xlim([0 800])
% legend('RMS','dx','dy','dz')
xlabel('Solution Number')
ylabel('Distance, [mm]')
% set(f1,'FontSize',24)
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
% set(gca,'FontSize',10,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',10,'fontWeight','bold')
hold off

%%



print(fig1, 'DatumCheck.eps','-depsc','-r300');
