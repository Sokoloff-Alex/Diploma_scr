%% Cross-section

Meridian = 13;
width = 0.5;
L = Meridian-width;
R = Meridian+width;
T = 49.0;
B = 45.5;
r = 1:198;

clear  LongG LatG V_prof V_profV

%% for Horizontal 

Outliers   = {'HELM', 'WIEN', 'FERR', 'FERH', 'OGAG', 'SOND', 'OBE2', ...
                'ROHR','BIWI','BI2I'};
%             
% Outliers   = {'HELM','WIEN','OGAG','OBE2', ...
%               'ROHR','BIWI','BI2I','MANS'};
            
%% for Vertical 
Outliers = {'ELMO','WIEN','FERR', ...
            'BIWI','BI2I','MANS','FFMJ','MOGN','WLBH', ...
            'TRF2','KRBG','OBE2','WT21','HKBL','PATK','PAT2', ...
            'HRIE','KTZ2', 'WLBH','ENTZ','OGAG'};
%%        
iiOut = selectRange(names, Outliers);
iiSel = setdiff(1:198, iiOut);

CovVenu2  = extractCovariance(CovVenu, iiSel, [1 2 3], 'no split');
V_enu_res = [Ve_res(iiSel), Vn_res(iiSel), Vu_res(iiSel)];
%% for 
CovScale = 20;
ii = (1:size(CovVenu2,1)/3)*3-2;
CovVenu3 = CovVenu2;
CovVenu3(ii,  ii)   = CovVenu3(ii,  ii)   * 10^2;
CovVenu3(ii+1,ii+1) = CovVenu3(ii+1,ii+1) * 10^2;
CovVenu3(ii+2,ii+2) = CovVenu3(ii+2,ii+2) * 10^2;


[LongG, LatG, V_profV,rmsFit_Vcs, V_SigPred_Vcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu3, ... 
    1, [Meridian Meridian], [B T], 0.125, 100, 5, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

[LongG, LatG, V_profH, rmsFit_Hcs, V_SigPred_Hcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu3, ... 
    1, [Meridian Meridian], [B T], 0.125, 100, 5, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

%%
iiB = r(lat >= B);
iiT = r(lat <= T);
iiR = r(long <= R);
iiL = r(long >= L);
iiBox = intersect(intersect(iiL, iiR), intersect(iiT,iiB) );
iiBox = setdiff(iiBox,iiOut);
iiBox = setdiff(iiBox, selectRange(names, {'NOVE', 'PORE','ALPE' } ));

CovVenuBox = extractCovariance(CovVenu, iiBox, [1 2 3], 'split');
stdVenu = sqrt(diag(CovVenuBox));

Vnerr = stdVenu((1:length(iiBox))*3-1)*1000 * 28;
Vuerr = stdVenu((1:length(iiBox))*3  )*1000 * 19;

Vnerr(Vnerr < 0.1) = 0.1;
Vuerr(Vuerr < 0.1) = 0.1;

%% topo
clc
nLon = size(ETOPO_Alps.Etopo_Europe,1);
nLat = size(ETOPO_Alps.Etopo_Europe,2);
step = ETOPO_Alps.refvec_Etopo(1);
Lat0 = ETOPO_Alps.refvec_Etopo(2);
Lon0 = ETOPO_Alps.refvec_Etopo(3);
clear elevProfile vLat elevProfile
elevProfile = ETOPO_Alps.Etopo_Europe(:, ((Meridian-width-Lon0)*step):((Meridian+width-Lon0)*step));
elevProfMean = mean(elevProfile,2);
elevProfMA   = movaverage(elevProfMean,5);
elevProfMax  = max(elevProfile')';
elevProfMin  = min(elevProfile')';
vLat = ( Lat0-(size(ETOPO_Alps.Etopo_Europe,1)-1)/step ):1/step:Lat0;


%% improve plot quality
%
% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;
clc

x_lim = [45 49];

% Defaults for this blog post
width = 4;    % Width in inches
height= 4;    % Height in inches
alw = 1;      % AxesLineWidth
fsz = 8;      % Fontsize
lw =  1;      % LineWidth
msz = 8;      % MarkerSize

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

V_SigPred_Vcs2 = abs(V_SigPred_Vcs)+0.05;
V_SigPred_Hcs2 = abs(V_SigPred_Hcs)+0.05;

close all
fig = figure(5);
hs(1) = subplot(9,1,1:3);
hold on
xv = [LatG;flipud(LatG)];
yv = [V_profH(:,3)*1000+V_SigPred_Vcs2(:,3); flipud(V_profH(:,3)*1000-V_SigPred_Vcs2(:,3))];
patch(xv, yv, 1, 'FaceColor', [240 0 0 ]/255 ,'EdgeColor','none' , 'FaceAlpha',0.2);
% pl1 = plot(LatG,V_profV(:,3)*1000+V_SigPred_Vcs(:,3),'--r','LineWidth',1);
% pl2 = plot(LatG,V_profV(:,3)*1000-V_SigPred_Vcs(:,3),'--r','LineWidth',1);
pl3 = plot(LatG,V_profV(:,3)*1000, 'LineWidth',2, 'Color',[170 0 0 ]/255);
pl4 = plot(lat(iiBox),Vu_res(iiBox)*1000,'ob', 'MarkerFaceColor','b', 'MarkerSize',3);
% text(lat(iiBox)+0.05,Vu_res(iiBox)*1000,names(iiBox))
errorbar(lat(iiBox),Vu_res(iiBox)*1000,Vuerr,'.b','LineWidth',.5);
% legend([pl3 pl1 pl4],'V_U LSC','LSC error','V_U obs')
legend([pl3 pl4],'V_U LSC','V_U obs')
title(['Velocity and topography profiles, along ',num2str(Meridian),'^oE Meridian'])
ylabel('V_U [mm/yr]')
xlim(x_lim)
ylim([-1.5 3])
set(gca, 'XTickLabel', [])
hold off

hs(2) = subplot(9,1,4:6);
hold on
xn = [LatG;flipud(LatG)];
yn = [V_profH(:,2)*1000+V_SigPred_Hcs2(:,2); flipud(V_profH(:,2)*1000-V_SigPred_Hcs2(:,2))];
patch(xn, yn, 1, 'FaceColor', [240 0 0 ]/255 ,'EdgeColor','none' , 'FaceAlpha',0.2);
% pl4 = plot(LatG,V_prof(:,2)*1000+V_SigPred_Hcs(:,2),'--r','LineWidth',1);
% pl5 = plot(LatG,V_prof(:,2)*1000-V_SigPred_Hcs(:,2),'--r','LineWidth',1);
pl6 = plot(LatG,V_profH(:,2)*1000, 'LineWidth',2, 'Color',[170 0 0 ]/255);
pl7 = plot(lat(iiBox),Vn_res(iiBox)*1000,'ob', 'MarkerFaceColor','b', 'MarkerSize',3);
% text(lat(iiBox)+0.05,Vn_res(iiBox)*1000,names(iiBox))
pl8 = errorbar(lat(iiBox),Vn_res(iiBox)*1000,Vnerr,'.b','LineWidth',.5);
% legend([pl6 pl4 pl7],'V_N LSC','LSC error','V_N obs')
legend([pl6 pl7],'V_N LSC','V_N obs')
xlim(x_lim)
ylim([-1 3])
set(gca, 'XTickLabel', [])
%title('North component')
ylabel('V_N [mm/yr]')
hold off

hs(3) = subplot(9,1,7:9);
hold on
pl2 = plot(vLat, elevProfMax/1000,'LineWidth',1,'Color',[.5 .5 .5]);
pl3 = plot(vLat, elevProfMin/1000,'LineWidth',1,'Color',[.5 .5 .5]);
pl1 = plot(vLat, elevProfMA/1000, 'LineWidth',3, 'Color',[34 139 34]/255);
x = [vLat, 50 24];
y1 = [elevProfMin; -3000; -3000 ]'/1000;
y2 = [elevProfMax; -3000; -3000 ]'/1000;

pl10 = patch(x, y2, 1, 'FaceColor', [.5 .5 .5]/2 ,'EdgeColor','none' , 'FaceAlpha',0.2);
pl11 = patch(x, y1, 1, 'FaceColor', [.5 .5 .5]/2 ,'EdgeColor','none' , 'FaceAlpha',0.2);
plot(x_lim, [0 0],'b')

legend([pl1 pl2 pl3 pl6 pl4 pl7],'Mean topography','topoMax','tomoMin')
% legend([pl1],'Mean topography')

xlim(x_lim)
ylim([-2.6, max(elevProfMax/1000)])
ylim([-2.6, max(elevProfMax/1000)])

xlabel('Latitude, [deg]')
ylabel('Topography [km]')
hold off

% close all


%% 
print(fig, '../dat/Pics/Cross-section_6.5meridian_4.eps','-depsc','-r300');
print(fig, '../dat/Pics/Cross-section_6.5meridian_4.pdf','-dpdf', '-r300');

% system('gs -o -q -sDEVICE=png256 -dEPSCrop -r300 -oimprovedExample_eps.png improvedExample.eps');

%%
clc
system('convert improvedExample.eps improvedExample.png');

%%  dave for GMT
GMT_dir = '~/Alpen_Check/MAP/Profiles/dat_13E/';

write_xyzTable([vLat', elevProfMax],  [GMT_dir,'topoMax.txt'],  '%5.2f %5.1f \n');
write_xyzTable([vLat', elevProfMin],  [GMT_dir,'topoMin.txt'],  '%5.2f %5.1f \n');
write_xyzTable([vLat', elevProfMean], [GMT_dir,'topoMean.txt'], '%5.2f %5.1f \n');

write_xyzTable([lat(iiBox),Vu_res(iiBox)*1000,Vuerr], [GMT_dir,'Vu_obs.txt'], '%7.4f %7.3f  %7.3f \n');
write_xyzTable([lat(iiBox),Vn_res(iiBox)*1000,Vnerr], [GMT_dir,'Vn_obs.txt'], '%7.4f %7.3f  %7.3f \n');

write_xyzTable([LatG,V_profV(:,3)*1000, V_SigPred_Vcs2(:,3)], [GMT_dir,'Vu_LSC.txt'], '%7.4f %7.3f  %7.3f \n');
write_xyzTable([LatG,V_profH(:,2)*1000, V_SigPred_Hcs2(:,2)], [GMT_dir,'Vn_LSC.txt'], '%7.4f %7.3f  %7.3f \n');

write_xyzTable([xn yn], [GMT_dir,'Vn_LSC_err.txt'], '%7.4f %7.3f \n');
write_xyzTable([xv yv], [GMT_dir,'Vu_LSC_err.txt'], '%7.4f %7.3f \n');




