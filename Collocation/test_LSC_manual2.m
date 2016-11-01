% test_LSC

close all
clear all
clc
%%

long_test = [-1 0 0 2 2    6 6 8 8  9]';
lat_test  = [ 1 0 2 0 2    0 2 0 2  1]';
% Vn_test   = [ 1 2 2 2 2    0 0 0 0 -1]';
% Ve_test   = [-1 0 0 0 0    2 2 2 2  1]';

Vn_test = [ 1.0054    1.9918    1.9893    1.9893    2.0087    0.0098   -0.0176   -0.0015    0.0025   -1.0061]';
Ve_test = [-1.0000    0.0033    0.0120   -0.0140    0.0021    1.9801    2.0003    2.0193    2.0094    0.9999]';
% Vn_test = Vn_test + randn(size(Vn_test))/100;
% Ve_test = Ve_test + randn(size(Ve_test))/100;

range = 1:length(lat_test);
V_obs = [Vn_test; Ve_test];
tic
lim = 3;
Max_Dist = 700;
range = 1:length(long_test);
p = 0;
clear V_pred_2 V_pred_3 V_pred_4 V_pred_5 LongGrid_3 LatGrid_3
for iLong = -2:0.5:12
    for iLat = -1:0.5:2
        arc = distance(iLat, iLong, lat_test, long_test) * 111 ; % km
        sel = range(arc < Max_Dist);  
        if length(sel) < lim % add more stations
            add = sort(arc);
            add = add(1:lim); % ad 2 sites, as 1st prop already included
            iadd = ismember(arc,add);
            inew = range(iadd);
            sel = unique(sort([sel, inew]));
        end
            p = p + 1;
            LongGrid_3(p) = iLong;
            LatGrid_3(p)  = iLat;   
            V_pred_3(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'exp1', '', '', 'bias', 'tail');
            V_pred_4(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'exp1', 'no corr','bias', 'tail');
                     
            % alternative
            V_pred_5(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'Hirvonen', '-v','','no corr');
            V_pred_6(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'Hirvonen' );
%           V_pred_6(:,p) = solve_LSC(iLat, iLong,lat_test, long_test, Vn_test/1000, Ve_test/1000,'Hirvonen', 'no plot')*1000;
    end
end
toc
disp('done')

%%
s = 0.2
clc
clr = lines(8);
try
  close (fig3)
end
fig3 = figure(3);
hold on; grid on
axis equal
xlim([-5 10])
ylim([-3 5])
plot(long_test,lat_test,'ok')
plot(LongGrid_3, LatGrid_3,'ok')
quiver(long_test,lat_test, Ve_test*s, Vn_test*s, 0, 'r','LineWidth',2)
% quiver(LongGrid_3, LatGrid_3, V_pred_2(2,:)*s, V_pred_2(1,:)*s, 0, 'Color', clr(1,:))
quiver(LongGrid_3, LatGrid_3, V_pred_3(2,:)*s, V_pred_3(1,:)*s, 0, 'b')
quiver(LongGrid_3, LatGrid_3, V_pred_4(2,:)*s, V_pred_4(1,:)*s, 0, 'g')
quiver(LongGrid_3, LatGrid_3, V_pred_5(2,:)*s, V_pred_5(1,:)*s, 0, 'm')
quiver(LongGrid_3, LatGrid_3, V_pred_6(2,:)*s, V_pred_6(1,:)*s, 0, 'k')

% for i = length(LatGrid_3)
%     ellipce_2D([km2deg(Max_Dist,6378*cosd(LatGrid_3(i))), Max_Dist/111], 0, [LongGrid_3(i), LatGrid_3(i)], 1)
% end

plot(long_test(sel), lat_test(sel), '*k')

 
 