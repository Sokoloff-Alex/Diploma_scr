function [fig] = plot_cross_section(point1, point2, width, ETOPO, Vobs, Vobs_err)
% function to plot cross-section between points (1 and 2)
% of toporgathy, and velocity field, with observations and
% interpolation (LSC) curve
% points specify by coordinates (lat, long) [deg, deg]
% width total in  [km]
%
% Alexandr Sokolov, KEG
% 11.01.2017

%% prepare ROI


%% Plot cross-sections

close all
fig = figure
subplot(2,1,1)
hold on
grid on
plot(vLat, elevProfMA/1000,'g','LineWidth',2)
plot(vLat, elevProfMax/1000,'--g','LineWidth',1)
plot(vLat, elevProfMin/1000,'--g','LineWidth',1)
plot(LatG,V_profV(:,3)*1000,'r','LineWidth',2)
plot(LatG,V_profV(:,3)*1000+V_SigPred_Vcs(:,3),'--r')
plot(LatG,V_profV(:,3)*1000-V_SigPred_Vcs(:,3),'--r')
plot(lat(iiBox),Vu_res(iiBox)*1000,'.b','LineWidth',1)
text(lat(iiBox)+0.05,Vu_res(iiBox)*1000,names(iiBox))
errorbar(lat(iiBox),Vu_res(iiBox)*1000,Vuerr,'.b')
% plot(vLat, elevProfMean/1000,'g','LineWidth',2)
legend('Topography Mean','topoMax','tomoMin','V_U LSC','LSE error','LSE error','V_U obs','ErrorBar 1 sigma')
title(['Vertical component, along ',num2str(Meridian),'^o Meridian'])
xlabel('Latitude, [deg]')
xlabel('Latitude, [deg]')
ylabel('V_U [mm/yr], Topo [km]')
xlim([B T])
subplot(2,1,2)
hold on
plot(vLat, elevProfMA/1000,'g','LineWidth',2)
plot(vLat, elevProfMax/1000,'--g','LineWidth',1)
plot(vLat, elevProfMin/1000,'--g','LineWidth',1)
plot(LatG,V_prof(:,2)*1000, 'r','LineWidth',2)
plot(LatG,V_prof(:,2)*1000+V_SigPred_Hcs(:,2),'--r')
plot(LatG,V_prof(:,2)*1000-V_SigPred_Hcs(:,2),'--r')
plot(lat(iiBox),Vn_res(iiBox)*1000,'.b')
text(lat(iiBox)+0.05,Vn_res(iiBox)*1000,names(iiBox))
errorbar(lat(iiBox),Vn_res(iiBox)*1000,Vnerr,'.b')
legend('Topo Mean','topoMax','tomoMin','V_N LSC','LSE error','LSE error','V_N obs','ErrorBar 1 sigma')
xlim([B T])
title('North component')
xlabel('Latitude, [deg]')
ylabel('V_N [mm/yr], Topo [km]')
grid on


end