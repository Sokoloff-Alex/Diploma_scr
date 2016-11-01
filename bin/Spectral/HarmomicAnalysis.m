function[c] = HarmomicAnalysis(values, time, TrendName)
% by Alexandr Sokolov

values = values*1000; % m -> mm
StatringEpoch = '2004-01-01';
time = time + datenum(StatringEpoch) - 1;

Year = 365.25; % days
% Year = 350; % days    % Also try Draconitic year = 350 days !!! 
MoonPeriod =  27.32;

T1 = Year / 1;   
T2 = Year / 2;   
T3 = Year / 3;   
T4 = Year / 4;   
T5 = Year / 6;


% T = [2 1 1/2 1/4 1/6]*Year/CyclePeriod;
% NumberOfHarmonics = length(T);
% 
% % make Desing Matrix
% TrendTerm = [ones(size(time,1),1),  time];
% HarmonicsTerms = zeros(size(time,1),length(T)*2);
% 
% for index = 1:2:NumberOfHarmonics*2
%     HarmonicsTerms(:,index)   = sin(2*pi/T((index+1)/2)*(time))';
%     HarmonicsTerms(:,index+1) = cos(2*pi/T((index+1)/2)*(time))';
% end
% 
% A = [TrendTerm HarmonicsTerms]; 


A = [ones(size(time,1),1)'; 
     time';
     sin(2*pi/T1*(time))';
     cos(2*pi/T1*(time))';
     sin(2*pi/T2*(time))'; 
     cos(2*pi/T2*(time))'; 
     sin(2*pi/T3*(time))'; 
     cos(2*pi/T3*(time))'; 
     sin(2*pi/T4*(time))'; 
     cos(2*pi/T4*(time))';
     sin(2*pi/T5*(time))';
     cos(2*pi/T5*(time))']';
 
c = (A'*A)^-1*A'*values;

error = A*c-values;
M = size(time,1);
N = 2;
variance = (error'*error)/(M-N);


y = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time)) + c(5)*sin(2*pi/T2*(time)) + c(6)*cos(2*pi/T2*(time)) + c(7)*sin(2*pi/T3*(time)) + c(8)*cos(2*pi/T3*(time)) + c(9)*sin(2*pi/T4*(time)) + c(10)*cos(2*pi/T4*(time)) + c(11)*sin(2*pi/T5*(time)) + c(12)*cos(2*pi/T5*(time));

TrendLine = c(1)+c(2)*(time);
f2 = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time));
f3 = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time)) + c(5)*sin(2*pi/T2*(time)) + c(6)*cos(2*pi/T2*(time));
f4 = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time)) + c(5)*sin(2*pi/T2*(time)) + c(6)*cos(2*pi/T2*(time)) + c(7)*sin(2*pi/T3*(time)) + c(8)*cos(2*pi/T3*(time));
f5 = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time)) + c(5)*sin(2*pi/T2*(time)) + c(6)*cos(2*pi/T2*(time)) + c(7)*sin(2*pi/T3*(time)) + c(8)*cos(2*pi/T3*(time))+ c(9)*sin(2*pi/T4*(time)) + c(10)*cos(2*pi/T4*(time));
f6 = c(1)+c(2)*(time) + c(3)*sin(2*pi/T1*(time)) + c(4)*cos(2*pi/T1*(time)) + c(5)*sin(2*pi/T2*(time)) + c(6)*cos(2*pi/T2*(time)) + c(7)*sin(2*pi/T3*(time)) + c(8)*cos(2*pi/T3*(time))+ c(9)*sin(2*pi/T4*(time)) + c(10)*cos(2*pi/T4*(time)) + + c(11)*sin(2*pi/T5*(time)) + c(12)*cos(2*pi/T5*(time));


h1 = c(3) * sin(2*pi/T1*(time)) + c(4) * cos(2*pi/T1*(time));
h2 = c(5) * sin(2*pi/T2*(time)) + c(6) * cos(2*pi/T2*(time));
h3 = c(7) * sin(2*pi/T3*(time)) + c(8) * cos(2*pi/T3*(time));
h4 = c(9) * sin(2*pi/T4*(time)) + c(10)* cos(2*pi/T4*(time));
h5 = c(11)* sin(2*pi/T5*(time)) + c(12)* cos(2*pi/T5*(time));

TimeSeries_no_harmonics = values - h1 - h2 - h3 - h4 - h5;

a_mmPerYear = c(2) * Year; %-> [mm/year]
SlopeTxt = ([' Resdual Trend: ',num2str(a_mmPerYear),' [mm/year]']);

if (~strcmp(TrendName,'no'))
    %%
    fig1 = figure;
    set(gcf,'PaperPositionMode','auto')
    set(fig1, 'Position', [0 0 1900 1000])
    subplot(4,2,1)
    hold on
    grid on
    plt_ts = plot(time,values,'.b')
    plt_tr =plot(time,TrendLine,'r','LineWidth',2)
    plt_h1= plot(time,f6,'m','LineWidth',2)
    legend(['timeseries, RMS = ', num2str(rms(values)), '; Mean = ',num2str(mean(values))],'trend','approximated')
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
    datetick('x','mmm-YYYY')
    
    subplot(4,2,2)
    hold on; grid on
    clr = colormap(lines(5))
    plot(time,h1,'LineWidth',2, 'Color',clr(1,:))
    plot(time,h2,'LineWidth',2, 'Color',clr(2,:))
    plot(time,h3,'LineWidth',2, 'Color',clr(3,:))
    plot(time,h4,'LineWidth',2, 'Color',clr(4,:))
    plot(time,h5,'LineWidth',2, 'Color',clr(5,:))
    legend(['1   year, Ampl = ',num2str(max(h1)),' [mm]'],...
           ['1/2 year, Ampl = ',num2str(max(h2)), '[mm]'],...
           ['1/3 year, Ampl = ',num2str(max(h3)), '[mm]'],...
           ['1/4 year, Ampl = ',num2str(max(h4)), '[mm]'],...
           ['1/6 year, Ampl = ',num2str(max(h5)), '[mm]'])
       
    ylabel('[mm]')
%     xlim([min(time) max(time)])
  datetick('x','mmm-YYYY')
    
    subplot(4,2,3)
    hold on; grid on
    plot(time,values-TrendLine,'.b')
    plot(time,h1,'r','LineWidth',2)
    legend(['timeseries minus trend, RMS = ', num2str(rms(values-TrendLine)) ],['annual, Amplidude = ',num2str(max(h1)), ' [mm]'])
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    
    subplot(4,2,4)
    hold on; grid on
    plot(time,values-f2,'.b')
    plot(time,h2,'r','LineWidth',2)
    legend(['timeseries - (trend + annual), RMS = ', num2str(rms(values-f2))],['1/2 year, Ampl = ',num2str(max(h2)), '[mm]'])
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    
    subplot(4,2,5)
    hold on, grid on
    plot(time,values-f3,'.b')
    plot(time,h3,'r','LineWidth',2)
    legend(['timeseries - (trend + annual + 1/2 year), RMS = ', num2str(rms(values-f3))],['1/3 year, Ampl = ',num2str(max(h3)), '[mm]'])
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    
    subplot(4,2,6)
    hold on; grid on
    plot(time,values-f4,'.b')
    plot(time,h4,'r','LineWidth',2)
    legend(['timeseries - (trend + annual + 1/2 + 1/3), RMS = ', num2str(rms(values-f4))],['1/4 year, Ampl = ',num2str(max(h4)), '[mm]'])
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    
    subplot(4,2,7)
    hold on; grid on
    plot(time,values-f5,'.b')
    plot(time,h5,'r','LineWidth',2)
    legend(['timeseries - (trend + annual + 1/2 + 1/3 + 1/4), RMS = ', num2str(rms(values-f5))],['1/6 year, Ampl = ',num2str(max(h5)), '[mm]'])
    xlabel('Time, [Year]')
    ylabel('[mm]')
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    subplot(4,2,8) 
    
    hold on; grid on
    plot(time,values -f6,'.b')
    plot(time,TrendLine,'r','LineWidth',2)
    xlabel('Time, [Year]')
    ylabel('[mm]')
    legend(['Timeseries without harmonics, RMS = ', num2str(rms(values -f6))], 'Trend')
    text('Position',[time(2), max(TimeSeries_no_harmonics)],'String',SlopeTxt)
    text('Position',[time(2), min(TimeSeries_no_harmonics)],'String',TrendName)
%     xlim([min(time) max(time)])
    ylim([-10 10])
  datetick('x','mmm-YYYY')
    % ylim([min(Trend_harmonics) max(Trend_harmonics)])

    
    %%
    fig2= figure;
    set(gcf,'PaperPositionMode','auto')
    set(fig2, 'Position', [0 0 1900 1000])
    subplot(2,1,1)
    hold on
    grid on
    plt_ts = plot(time,values,'.b')
    plt_h = plot(time,f6,'r','LineWidth',2)
    xlabel('time')
    ylabel('[mm]')
    legend(['timeseries, RMS = ', num2str(rms(values))],'approximation')
    datetick('x','mmm-YYYY')
    
    subplot(2,1,2)
    hold on
    grid on
    plt_1 = plot(time,TimeSeries_no_harmonics,'.b')
    plt_2 = plot(time,TrendLine,'r','LineWidth',2)
    xlabel('time')
    ylabel('[mm]')
    legend(['Residusals, RMS =', num2str(rms(TimeSeries_no_harmonics))], 'Trend')      
    datetick('x','mmm-YYYY')
    text('Position',[time(2), max(TimeSeries_no_harmonics)],'String',SlopeTxt)
    text('Position',[time(2), min(TimeSeries_no_harmonics)],'String',TrendName) 
    
%     print(fig1, '-dpng',['Results\Trends\Trend_All_',TrendName,'.png']);
%     print(fig2, '-dpng',['Results\Trends\Trend_',TrendName,'.png']);
end

end