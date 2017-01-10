%% Cross-section
Meridian = 13;
L = Meridian-0.5;
R = Meridian+0.5;
T = 49.0;
B = 45.5;
r = 1:198;
iiB = r(lat >= B);
iiT = r(lat <= T);
iiR = r(long <= R);
iiL = r(long >= L);
iiBox = intersect(intersect(iiL, iiR), intersect(iiT,iiB) );
iiBox = setdiff(iiBox,iiOut);
iiBox = setdiff(iiBox, selectRange(names, {'NOVE', 'PORE','ALPE' } ));
clear  LongG LatG V_prof V_profV
[LongG, LatG, V_profV,rmsFit_Vcs, V_SigPred_Vcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [Meridian Meridian], [B T], 0.125, 150, 5,  'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');
[LongG, LatG, V_prof, rmsFit_Hcs, V_SigPred_Hcs] = run_Collocation(long(iiSel), lat(iiSel), V_enu_res, CovVenu2, [Meridian Meridian], [B T], 0.125, 100, 5,  'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');

CovVenuBox = extractCovariance(CovVenu, iiBox, [1 2 3], 'split');
CavVenu = sqrt(diag(CovVenuBox));
Vuerr = CavVenu((1:length(iiBox))*3  )*1000*sqrt(1.8)*20;
Vnerr = CavVenu((1:length(iiBox))*3-2)*1000*sqrt(1.8)*20;
Vuerr(Vuerr < 0.1) = 0.1;
Vnerr(Vnerr < 0.1) = 0.1;

%%
clc
nLon = size(ETOPO_Alps.Etopo_Europe,1);
nLat = size(ETOPO_Alps.Etopo_Europe,2);
step = ETOPO_Alps.refvec_Etopo(1);
Lat0 = ETOPO_Alps.refvec_Etopo(2);
Lon0 = ETOPO_Alps.refvec_Etopo(3);
LonRange = [ Lon0, Lon0 + nLat/step];
LatRange = [ Lat0, Lat0 - nLon/step];
clear elevProfile vLat
elevProfile = ETOPO_Alps.Etopo_Europe(:, ((Meridian-0.5-Lon0)*step):((Meridian+0.5-Lon0)*step));
elevProfile = mean(elevProfile,2);
vLat = ( Lat0-(size(ETOPO_Alps.Etopo_Europe,1)-1)/step ):1/step:Lat0;

%
close all
figure(5)
subplot(2,1,1)
hold on
grid on
plot(vLat, elevProfile/1000,'g','LineWidth',2)
plot(LatG,V_profV(:,3)*1000,'r','LineWidth',2)
plot(LatG,V_profV(:,3)*1000+V_SigPred_Vcs(:,3),'--r')
plot(LatG,V_profV(:,3)*1000-V_SigPred_Vcs(:,3),'--r')
plot(lat(iiBox),Vu_res(iiBox)*1000,'.b','LineWidth',1)
text(lat(iiBox)+0.05,Vu_res(iiBox)*1000,names(iiBox))
errorbar(lat(iiBox),Vu_res(iiBox)*1000,Vuerr,'.b')
legend('Topography','V_U LSC','LSE error','LSE error','V_U obs')
title(['Vertical component, along ',num2str(Meridian),'^o Meridian'])
xlabel('Latitude, [deg]')
xlabel('Latitude, [deg]')
ylabel('V_U [mm/yr], Topo [km]')
xlim([B T])
subplot(2,1,2)
hold on
plot(vLat, elevProfile/1000,'g','LineWidth',2)
plot(LatG,V_prof(:,2)*1000, 'r','LineWidth',2)
plot(LatG,V_prof(:,2)*1000+V_SigPred_Hcs(:,2),'--r')
plot(LatG,V_prof(:,2)*1000-V_SigPred_Hcs(:,2),'--r')

plot(lat(iiBox),Vn_res(iiBox)*1000,'.b')
text(lat(iiBox)+0.05,Vn_res(iiBox)*1000,names(iiBox))
errorbar(lat(iiBox),Vn_res(iiBox)*1000,Vnerr,'.b')
xlim([B T])
title('North component')
legend('Topography','V_N LSC','LSE error','LSE error','V_N obs')
xlabel('Latitude, [deg]')
ylabel('V_N [mm/yr], Topo [km]')
grid on

