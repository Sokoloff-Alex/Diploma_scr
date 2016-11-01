
close all
clc

[Stations, Radoms, Records] = parse_OUT_file('/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FCS_out4.OUT');

%%
write_residuals_files(Records);
%%
% write_LLH(Records);
%%
 
read_PLT('/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FCS_out4.PLT', Stations)

%%
filename = '/home/gast/GPSDATA/CAMPAIGN52/ALP_NET/OUT/FCS_out4.PLT'
Station = {'ESAB'};
[MJD, res_N_stack, res_E_stack, res_U_stack, res_H_stack, res_3D_stack, epochNmr] = get_PLT_residuals(filename, Station);

%%
close all; clc
[c] = HarmomicAnalysis(res_E_stack, epochNmr, Station)

%%
close all; clc
fftAnalysis(epochNmr, res_E_stack*1000, 1/500, 1/340)
%%
figure(2)
periodogram(res_E_stack*1000, 'power')
ax = gca;
set(ax, 'XTick',[ 1/Year, 1/Year*2, 1/Year*3, 1/Year*4, 1/Year*6,1/Year*8, 1/Year*12]);
set(ax,'XTickLabel', {'1','1/2','1/3','1/4','1/6','1/8','1/12'});
xlabel(['Period, [years]; 1 year = ',num2str(Year),' days']);
xlim([0 1/10]);

%%
close all
figure(1)
subplot(3,1,1)
hold on
grid on
title(['Residuals, ',Station])
plot(epochNmr, res_N_stack*1000,'.b')
legend('North')
ylabel('[mm]')
ylim([-15 15])
subplot(3,1,2)
hold on
grid on
plot(epochNmr, res_E_stack*1000,'.b')
legend('East')
ylabel('[mm]')
ylim([-15 15])
subplot(3,1,3)
hold on
grid on
plot(epochNmr,res_U_stack*1000,'.b')
legend('Up')
ylabel('[mm]')
ylim([-30 30])



