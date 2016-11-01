% try FFT and harmonic analisys

%%
close all
% days = cell2mat(North(:,1));
% Res= cell2mat(North(:,2));
% days = cell2mat(East(:,1));
% Res= cell2mat(East(:,2));
days = cell2mat(q3(:,1));
Res= cell2mat(q3(:,2));
% days = cell2mat(Up(:,1));
% Res= cell2mat(Up(:,2));

%%
close all
figure(1)
hold on
grid on
plot(days-days(1), Res, '.b')
datetick('x','YYYY') 
ylabel('Residual North, [mm]')

%%
close all
fftAnalysis(days - days(1),Res, 0, 1/180)
%%

[f,P,prob] = lomb(Res,days - days(1) ,4,1/10);
%%
close all
figure(2)
hold on
grid on
   plot(f,P)
   [Pmax,jmax] = max(P)
   disp(['Most significant period is ',num2str(1/f(jmax)),...
        ' with FAP of ',num2str(prob(jmax))])

%%
close all
HarmomicAnalysis(Res,days, 'WELS, North')