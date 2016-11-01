function fftAnalysis(time, values, fmin, fmax)
%% FFT transform and filtering


StatringEpoch = '2004-01-01';

% replace missing with 0
counter = 1;
for i = 1:(time(end) - time(1))
    time_new(i) = i + time(1) - 1;
   if time(counter) == i+time(1)-1;
       counter = counter + 1;       
       values_new(i) = values(counter);
   else
        values_new(i) = 0;
%       values_new(i) = mean([values(counter-1:counter)]);
   end
end

%%

Year = 365.25; % days
% Year = 350 days; try draconitic year
MoonPeriod = 27.321662; % [days], Moon Sidereal month

% Sampling
Period=1; % day
SamplingFrequency = 1/Period; % 1/day
Duration = time(end) - time(1)*1/SamplingFrequency;
% time=0 : Period : Duration;            time = time(1:end-1);
frequency = -SamplingFrequency/2 : 1/Duration : SamplingFrequency/2; 
frequency = frequency(1:end-1);


Signal = values_new - nanmean(values_new);

% FFT
SIG = fft(Signal); SIG_sort = abs(fftshift(SIG));

% Filter
FIL_LP_out = rectpuls(frequency, fmax*2);
FIL_LP_in  = rectpuls(frequency, fmin*2);
FIL_sort = FIL_LP_out - FIL_LP_in;

FIL = ifftshift(FIL_sort); 

fil = ifft(FIL);
fil_sort = fftshift(fil);

% Filtering
SIG2 = SIG.*FIL;
SIG2_sort = abs(fftshift(SIG2));

% iFFT
SignalFiltered = ifft(SIG2);

time_new = time_new + datenum(StatringEpoch) - 1;
%% plots,   in time domain
figure
subplot(2,2,1)
hold on
plot(time_new, Signal,'.b');
plot(time_new,SignalFiltered,'r','LineWidth',2)
grid on
title('Time domain');
xlabel('time [yyyy]');
ylabel('Amplitude');
datetick('x','YYYY')
% ylim([-15 15])
legend(['RMS = ', num2str(rms(Signal)), ' ; Mean = ', num2str(mean(values)) ])

subplot(2,2,3)
plot(time_new,Signal - SignalFiltered,'.b')
grid on;
title('Signal - Filtered');
ylabel('Amplitude');
xlabel('time [yyyy]');
datetick('x','YYYY')
% ylim([-10 10])
legend(['RMS = ', num2str(rms(Signal - SignalFiltered))])


% in freq domain
subplot(2,2,2)
hold on; grid on
plot(frequency, SIG_sort);
plot(frequency, FIL_sort*max(SIG_sort)*1.1,'r','LineWidth',2);
% xlim([0 SamplingFrequency/2]);
xlim([0 1/30]);
title('Frequency  domain');
ylabel('Amplitude');

plot([1/(1/1*Year),  1/(1/1*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/2*Year),  1/(1/2*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/3*Year),  1/(1/3*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/4*Year),  1/(1/4*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/6*Year),  1/(1/6*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/8*Year),  1/(1/8*Year)],  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot([1/(1/12*Year), 1/(1/12*Year)], [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
% text('Position',[1/Year, max(SIG_sort)],'String','Annual')
ax = gca;
set(ax, 'XTick',[0, 1/Year, 1/Year*2, 1/Year*3, 1/Year*4, 1/Year*6,1/Year*8, 1/Year*12]);
set(ax,'XTickLabel', {'0','1','1/2','1/3','1/4','1/6','1/8','1/12'});
xlabel(['Period, [years]; 1 year = ',num2str(Year),' days']);
freq_pos  = frequency([ round(length(frequency)/2) : end]);
signal_sort_pos = SIG_sort([round(length(frequency)/2):end]);

[Max_Signal,Index] = max(signal_sort_pos); 
Period_Max = 1/abs(freq_pos(Index)); % Period of max peak, [days]
% text('Position',[1/Period_Max, max(SIG_sort)*0.9],'String',['Max Component = 1/',num2str(Period_Max)])
legend(['Max Component = 1/',num2str(Period_Max), ' day^-^1'])

hold off
% in freq domain
subplot(2,2,4)
hold on; grid on
freq_pos  = frequency([ round(length(frequency)/2) : end]);
signal_sort_pos = SIG_sort([round(length(frequency)/2):end]);
plot(log(freq_pos), signal_sort_pos);
plot(log(frequency), FIL_sort*max(signal_sort_pos)*1.1,'r','LineWidth',2);

plot(log([1/(2*Year),    1/(2*Year)]),    [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1*Year),    1/(1*Year)]),    [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/2*Year),  1/(1/2*Year)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/3*Year),  1/(1/3*Year)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/4*Year),  1/(1/4*Year)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/6*Year),  1/(1/6*Year)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/8*Year),  1/(1/8*Year)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
plot(log([1/(1/12*Year), 1/(1/12*Year)]), [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.5 0.5], 'LineWidth',2)
% plot(log([1/MoonPeriod,  1/MoonPeriod]),  [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.7 0.7], 'LineWidth',2)

plot(log([1/(1*Period_Max),    1/(1*Period_Max)]),    [0 max(SIG_sort)*1.05],'--', 'Color',[138,43,226]/256, 'LineWidth',2)
plot(log([1/(1/2*Period_Max),  1/(1/2*Period_Max)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[138,43,226]/256, 'LineWidth',2)
plot(log([1/(1/3*Period_Max),  1/(1/3*Period_Max)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[138,43,226]/256, 'LineWidth',2)
plot(log([1/(1/4*Period_Max),  1/(1/4*Period_Max)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[138,43,226]/256, 'LineWidth',2)
plot(log([1/(1/6*Period_Max),  1/(1/6*Period_Max)]),  [0 max(SIG_sort)*1.05],'--', 'Color',[138,43,226]/256, 'LineWidth',2)



% plot(log([1/7,           1/7]),           [0 max(SIG_sort)*1.05],'--', 'Color',[0 0.7 0.7], 'LineWidth',2)
ax = gca;
set(ax, 'XTick', log([1/Year*0.5, 1/Year, 1/Year*2, 1/Year*3, 1/Year*4, 1/Year*6,1/Year*8, 1/Year*12]));
set(ax,'XTickLabel', {'2', '1','1/2','1/3','1/4','1/6','1/8','1/12'});
xlim(log([1/(3*Year) 1/15]));
xlabel(['Period, [years]; 1 year = ',num2str(Year),' days']);
title('Frequency  domain');
ylabel('Amplitude');
% text('Position',[log(1/Period_Max), max(SIG_sort)*0.9],'String',['Max Component = 1/',num2str(Period_Max)])
legend(['Max Component = 1/',num2str(Period_Max), ' day^-^1'])
hold off

%%
end