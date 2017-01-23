%% Cross-section
Meridian = 6.5;
width = 0.5;
L = Meridian-width;
R = Meridian+width;
T = 49.0;
B = 42.75;
r = 1:198;
iiB = r(lat >= B);
iiT = r(lat <= T);
iiR = r(long <= R);
iiL = r(long >= L);
iiBox = intersect(intersect(iiL, iiR), intersect(iiT,iiB) );
iiBox = setdiff(iiBox,iiOut);
iiBox = setdiff(iiBox, selectRange(names, {'NOVE', 'PORE','ALPE' } ));
clear  LongG LatG V_prof V_profV
[LongG, LatG, V_profV,rmsFit_Vcs, V_SigPred_Vcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [Meridian Meridian], [B T], 0.125, 100, 5,  'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');
[LongG, LatG, V_prof, rmsFit_Hcs, V_SigPred_Hcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [Meridian Meridian], [B T], 0.125, 100, 5,  'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

CovVenuBox = extractCovariance(CovVenu, iiBox, [1 2 3], 'split');
CavVenu = sqrt(diag(CovVenuBox));
Vuerr = CavVenu((1:length(iiBox))*3  )*1000*sqrt(1.8)*20;
Vnerr = CavVenu((1:length(iiBox))*3-2)*1000*sqrt(1.8)*20;
Vuerr(Vuerr < 0.1) = 0.1;
Vnerr(Vnerr < 0.1) = 0.1;

%
clc
nLon = size(ETOPO_Alps.Etopo_Europe,1);
nLat = size(ETOPO_Alps.Etopo_Europe,2);
step = ETOPO_Alps.refvec_Etopo(1);
Lat0 = ETOPO_Alps.refvec_Etopo(2);
Lon0 = ETOPO_Alps.refvec_Etopo(3);
LonRange = [ Lon0, Lon0 + nLat/step];
LatRange = [ Lat0, Lat0 - nLon/step];
clear elevProfile vLat
elevProfile = ETOPO_Alps.Etopo_Europe(:, ((Meridian-width-Lon0)*step):((Meridian+width-Lon0)*step));
elevProfMean = mean(elevProfile,2);
elevProfMA   = movaverage(elevProfMean,5);
elevProfMax  = max( elevProfile')';
elevProfMin  = min( elevProfile')';
vLat = ( Lat0-(size(ETOPO_Alps.Etopo_Europe,1)-1)/step ):1/step:Lat0;


%% improve plot quality
%
% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;
clc

% Defaults for this blog post
width = 7;     % Width in inches
height = 7;    % Height in inches
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
fig = figure(5);
hs(1) = subplot(11,1,[1:4]);
hold on
grid on
% x = [LatG;flipud(LatG)];
% y = [V_prof(:,2)*1000+V_SigPred_Hcs(:,2); flipud(V_prof(:,2)*1000-V_SigPred_Hcs(:,2))];
% pl9 = patch(x, y, 1, 'FaceColor', [240 0 0 ]/255 ,'EdgeColor','none' , 'FaceAlpha',0.2)
pl4 = plot(LatG,V_profV(:,3)*1000+V_SigPred_Vcs(:,3),'--r','LineWidth',.25);
pl5 = plot(LatG,V_profV(:,3)*1000-V_SigPred_Vcs(:,3),'--r','LineWidth',.25);
pl6 = plot(LatG,V_profV(:,3)*1000, 'LineWidth',2, 'Color',[170 0 0 ]/255);
pl7 = plot(lat(iiBox),Vu_res(iiBox)*1000,'ob', 'MarkerFaceColor','b', 'MarkerSize',3);
% text(lat(iiBox)+0.05,Vu_res(iiBox)*1000,names(iiBox))
pl8 = errorbar(lat(iiBox),Vu_res(iiBox)*1000,Vuerr,'.b','LineWidth',.5);
legend([pl6 pl4 pl7],'V_U LSC','LSC error','V_U obs')
title(['Velocity and topography profiles, along ',num2str(Meridian),'^oE Meridian'])
ylabel('V_U [mm/yr]')
xlim([B T])
set(gca, 'XTickLabel', [])

hs(2) = subplot(11,1,[5:8]);
hold on
% x = [LatG;flipud(LatG)];
% y = [V_prof(:,2)*1000+V_SigPred_Hcs(:,2); flipud(V_prof(:,2)*1000-V_SigPred_Hcs(:,2))];
% pl9 = patch(x, y, 1, 'FaceColor', [240 0 0 ]/255 ,'EdgeColor','none' , 'FaceAlpha',0.2)
pl4 = plot(LatG,V_prof(:,2)*1000+V_SigPred_Hcs(:,2),'--r','LineWidth',.25);
pl5 = plot(LatG,V_prof(:,2)*1000-V_SigPred_Hcs(:,2),'--r','LineWidth',.25);
pl6 = plot(LatG,V_prof(:,2)*1000, 'LineWidth',2, 'Color',[170 0 0 ]/255);
pl7 = plot(lat(iiBox),Vn_res(iiBox)*1000,'ob', 'MarkerFaceColor','b', 'MarkerSize',3);
% text(lat(iiBox)+0.05,Vn_res(iiBox)*1000,names(iiBox))
pl8 = errorbar(lat(iiBox),Vn_res(iiBox)*1000,Vnerr,'.b','LineWidth',.5);
legend([pl6 pl4 pl7],'V_N LSC','LSC error','V_N obs')
xlim([B T])
ylim([-1 1])
set(gca, 'XTickLabel', [])
%title('North component')
ylabel('V_N [mm/yr]')
grid on
hold off

hs(3) = subplot(11,1,[9:11]);
hold on
pl2 = plot(vLat, elevProfMax/1000,'LineWidth',.5,'Color',[.5 .5 .5]);
pl3 = plot(vLat, elevProfMin/1000,'LineWidth',.5,'Color',[.5 .5 .5]);
pl1 = plot(vLat, elevProfMA/1000, 'LineWidth',2, 'Color',[34 139 34]/255);
legend([pl1 pl2 pl3 pl6 pl4 pl7],'Mean topography','topoMax','tomoMin')
xlim([B T])
ylim([-2.6, max(elevProfMax/1000)])
xlabel('Latitude, [deg]')
ylabel('Topography [km]')
grid on
hold off

% close all


%% 
print(fig, 'dat/Pics/Cross-section_6.5meridian_3.eps','-depsc','-r300');
% system('gs -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png improvedExample.eps');

%%
clc
system('convert improvedExample.eps improvedExample.png');





