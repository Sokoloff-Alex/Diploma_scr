% Oscillations

T = PERI365;

names_t = cell2mat(T.ACCE);
names_t = names_t(:,2:5);
names_all

r = [1:length(names_all)]';

ii = r(ismember(names_all, names_t))

%%
clc
close all


%%

ii = ii(1:274);

clc
close all
figure(2)
colormap('hsv')
subplot(2,2,1); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName8*10000+1, PERP365.VarName8 ,'filled' )
title('East')
subplot(2,2,2); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName9*10000+1, PERP365.VarName9 ,'filled' )
title('North')
subplot(2,2,3); hold on
scatter(wrapTo180(long_all(ii)), lat_all(ii), PERI365.VarName10*10000+1, PERP365.VarName10 ,'filled' )
title('Up')
