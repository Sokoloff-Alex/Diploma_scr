% test_LSC

close all
% clear all
clc

%%
clc
clear LongGrid LatGrid  Vn_pred_stack Ve_pred_stack Vu_pred_stack

range = 1:length(lat);
Max_Dist = 200;

tic
p = 0;
for iLong = 7:0.5:12
    for iLat = 47:0.5:49
        p = p + 1;
        LongGrid(p) = iLong;
        LatGrid(p)  = iLat;      

        arc = greatcircleArc(iLat, iLong,lat, long) * 111 ; % km
        sel = range(arc < Max_Dist);    
        sel = intersect(sel, Selected);
        
        V_pred = solve_LSC_3D(iLat, iLong, lat(sel), long(sel), Ve_res(sel), Vn_res(sel), Vu_res(sel), Max_Dist); % [mm/yr]
        
        Vn_pred_stack(p) = V_pred(1);
        Ve_pred_stack(p) = V_pred(2);
        Vu_pred_stack(p) = V_pred(3);
    end
end
t = toc


%%
clc
s = 1000;
close all
figure(3)
hold on
grid on
% plot(long,lat,'or')
plot(LongGrid, LatGrid,'ok')
% error_ellipse([Max_Dist/111, 0 ; 0 ,Max_Dist/111], [iLong, iLat])
% quiver(long(Selected),lat(Selected), Ve_res(Selected)*s, Vn_res(Selected)*s, 0, 'r')
quiver3(long(Selected),lat(Selected), zeros(length(Selected),1), Vu_res(Selected)*s, 0, 'g')
% quiver(LongGrid, LatGrid, Ve_pred_stack*s, Vn_pred_stack*s,0,'b')
quiver3(LongGrid, LatGrid, zeros(1,length(LongGrid)), Vu_pred_stack*s,0,'m')
xlim([min(LongGrid)-2 max(LongGrid)+2])
ylim([min(LatGrid)-2  max(LatGrid)+2])

%%
clc
close all
figure(3)
hold on ; grid on
% plot3(LongGrid, LatGrid,'ok')
% error_ellipse([Max_Dist/111, 0 ; 0 ,Max_Dist/111], [iLong, iLat])
% quiver(long(Selected),lat(Selected), Ve_res(Selected)*s, Vn_res(Selected)*s, 0, 'r')
quiver3(long(Selected),lat(Selected),zeros(length(Selected),1),Ve_res(Selected)*s, Vn_res(Selected)*s, Vu_res(Selected)*s, 0)
% quiver(LongGrid, LatGrid, Ve_pred_stack*s, Vn_pred_stack*s,0,'b')
quiver3(LongGrid, LatGrid, zeros(1,length(LongGrid)),Ve_pred_stack*s, Vn_pred_stack*s, Vu_pred_stack*s,0)
xlim([min(LongGrid)-2 max(LongGrid)+2])
ylim([min(LatGrid)-2  max(LatGrid)+2])

 
 
 