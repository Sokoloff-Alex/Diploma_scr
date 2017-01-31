% data_continuity

close all
clear all
clc


load('../dat/SNX/SINEX.mat')

%% compute common observation period
% 
t_start = SINEX.SOLUTION.EPOCHS.DATA_START;
t_end   = SINEX.SOLUTION.EPOCHS.DATA_END;

t_start = [str2num(t_start(:,1:2)) + str2num(t_start(:,4:6))/365.25]
t_end   = [str2num(t_end(  :,1:2)) + str2num(t_end(  :,4:6))/365.25]

dt = t_end - t_start;
% 
% %%
[dtobs_sum] = merge_stations_sum(dt,names_all);

%% plot histogramm

close all
figure(1)
hold on
grid on

[nelements1,centers1] = hist(dtobs_sum, [0:1:13])
[nelements2,centers2] = hist(dt,        [0:1:13])

bar(centers1, nelements1, 'b', 0.75)
bar(centers2, nelements2, 'r', 0.5)

legend('total time for each stations','Artificial','location','NorthWest')
xlabel('Duration, years')
ylabel('# of stations')

