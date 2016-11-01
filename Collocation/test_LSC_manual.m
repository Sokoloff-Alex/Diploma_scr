% test_LSC
close all
clear all
clc

%%

long_test = [0 0 2 4 5]';
lat_test  = [4 0 1 3 2]';
Vn_test   = [1 3 1 2 3]';
Ve_test   = [2 1 3 3 1]';

range = 1:length(lat_test);
V_obs = [Vn_test; Ve_test];
tic
lim = 3;
Max_Dist = 700;
range = 1:length(lat_test);
p = 0;
clear V_pred_2 V_pred_3 V_pred_4 V_pred_5 LongGrid_3 LatGrid_3
for iLong =  -5:10
    for iLat = -2:7
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
          % V_pred_2(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel)/1000, Ve_test(sel)/1000, 'with plot')*1000;
            V_pred_3(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'exp1',  'bias',  'tail');
            V_pred_4(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'exp1', 'no corr','bias', 'tail');
            V_pred_5(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'Hirvonen', '-v');
            V_pred_6(:,p) = solve_LSC(iLat, iLong, lat_test(sel), long_test(sel), Vn_test(sel), Ve_test(sel),'Hirvonen', 'no corr');
                                  
            % alternative
%           V_pred_5(:,p) = solve_LSC_alt(iLat, iLong,lat_test, long_test, Vn_test,      Ve_test,      'with plot');
%           V_pred_6(:,p) = solve_LSC_alt(iLat, iLong,lat_test, long_test, Vn_test/1000, Ve_test/1000, 'no plot')*1000;
    end
end
toc
disp('done')

%%
s = 0.25
clr = lines(8);
try
  close (fig3)
end
fig3 = figure(3);
hold on; grid on
axis equal
xlim([-7 12])
ylim([-5 8])
plot(long_test,lat_test,'ok')
plot(LongGrid_3, LatGrid_3,'ok')
quiver(long_test,lat_test, Ve_test*s, Vn_test*s, 0, 'r','LineWidth',2)
% quiver(LongGrid_3, LatGrid_3, V_pred_2(2,:)*s, V_pred_2(1,:)*s, 0, 'Color', clr(1,:))
quiver(LongGrid_3, LatGrid_3, V_pred_3(2,:)*s, V_pred_3(1,:)*s, 0, 'b')
quiver(LongGrid_3, LatGrid_3, V_pred_4(2,:)*s, V_pred_4(1,:)*s, 0, 'g')
% quiver(LongGrid_3, LatGrid_3, V_pred_5(2,:)*s, V_pred_5(1,:)*s, 0, 'm')
% quiver(LongGrid_3, LatGrid_3, V_pred_6(2,:)*s, V_pred_6(1,:)*s, 0, 'k')

ellipce_2D([km2deg(Max_Dist,6378*cosd(4)), Max_Dist/111], 0, [4, 5], 1)
plot(long_test(sel), lat_test(sel), '*k')

 
 