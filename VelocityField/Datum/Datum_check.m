% Datum check

Datum_dr = DatumcheckFMCIGBW7;

%%
clc
[XYZ] = lla2ecef([45, 8, 0]);
[de, dn, du] = XYZ2ENU(repmat(XYZ,length(Datum_dr),1), Datum_dr(:,3:5));

close all
f1 = figure(1);
hold on
grid on
% plot(Datum_dr(:,1),Datum_dr(:,2)*1000,'.k')
plot(Datum_dr(:,1),Datum_dr(:,3)*1000+10,'.r')
plot(Datum_dr(:,1),Datum_dr(:,4)*1000+20,'.g')
plot(Datum_dr(:,1),Datum_dr(:,5)*1000+30,'.b')
plot(Datum_dr(:,1), de*1000+40, '.r')
plot(Datum_dr(:,1), dn*1000+50, '.g')
plot(Datum_dr(:,1), du*1000+60, '.b')
% text(660,2,   'RMS')
text(660,1+10, ['dX, rms = ', num2str( rms(Datum_dr(:,3))*1000,'%3.2f')])
text(660,1+20, ['dY, rms = ', num2str( rms(Datum_dr(:,4))*1000,'%3.2f')])
text(660,1+30, ['dZ, rms = ', num2str( rms(Datum_dr(:,5))*1000,'%3.2f')])
text(660,1+40, ['dE, rms = ', num2str( rms(de)*1000,'%3.2f')])
text(660,1+50, ['dN, rms = ', num2str( rms(dn)*1000,'%3.2f')])
text(660,1+60, ['dU, rms = ', num2str( rms(du)*1000,'%3.2f')])
xlim([0 800])
% legend('RMS','dx','dy','dz')
xlabel('Solution Number')
ylabel('Distance, [mm]')
% set(f1,'FontSize',24)
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
set(gca,'FontSize',30,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',30,'fontWeight','bold')
hold off
