% SAPOS LSC

clc
iiSel=[1:length(long_all)];
V_enu_res = [Ve_res_all(iiSel), Vn_res_all(iiSel), Vu_res_all(iiSel)];

CovVenu2 = extractCovariance(CovVenuSNX, iiSel, [1 2 3], 'no split');
Cov_scale = 20;


%%
clc
[LongGrid, LatGrid, V_def_tr1, rmsFit, V_Sig_tr] = run_Collocation(long_all(iiSel), lat_all(iiSel), V_enu_res, CovVenuSNX, Cov_scale, [7 14], [46 51], 0.25, 80, 7, 'exp1', '-v', 'bias', 'tail 0', 'no corr', 'filter');



%%  Vertical
try
    close(fig1)
end
fig1 = figure(1);
hold on
sc=300
text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all, zeros(size(Vu_res_all)),      Vu_res_all*sc,    0, 'b');
quiver(LongGrid,            LatGrid, zeros(size(V_def_tr1(:,1))),  V_def_tr1(:,3)*sc,0, 'r', 'lineWidth',0.5)
xlim([6 15])
ylim([45 52])


%%  Horizontal
try
    close(fig2)
end
fig2 = figure(2);
hold on
sc=500
text(wrapTo180(long_all),   lat_all, names_all);
quiver(wrapTo180(long_all), lat_all, Ve_res_all*sc,      Vn_res_all*sc,    0, 'b');
quiver(LongGrid,            LatGrid, V_def_tr1(:,1)*sc,V_def_tr1(:,2)*sc,0, 'r', 'lineWidth',0.5)
xlim([6 15])
ylim([45 52])

