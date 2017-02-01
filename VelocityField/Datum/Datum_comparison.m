% Datum_comparison
%
% compare CRD and VEL between datums: 
% EPN EPN_A_IGb08.SNX
% IGS IGS16P26.SNX
% IGb IGb08.SNX
%
% for sations: GRAS, GRAZ, LROC, MEDI, WTZR, ZIMM

close all
clear all
clc
%%
% 
% ALP_SNX = readSNX('../../dat/SNX/FMC_IGB_W7_noCov.SNX','no COVA');
IGB_SNX = readSNX('../../dat/SNX/IGb08_noCov.snx', 'no COVA');
IGS_SNX = readSNX('../../dat/SNX/IGS16P26.snx', 'no COVA');
EPN_SNX = readSNX('../../dat/SNX/EPN_A_IGb08_no_COVA.SNX', 'no COVA');

% save('../../dat/SNX/IGb08_noCov.SNX.mat',  'IGB_SNX')
% save('../../dat/SNX/IGS16P26.SNX.mat',     'IGS_SNX')
% save('../../dat/SNX/IGb08_no_COVA.SNX.mat','EPN_SNX')

%%
load('../../dat/SNX/IGb08_noCov.SNX.mat')
load('../../dat/SNX/IGSIGS16P26.SNX.mat')
load('../../dat/SNX/IGb08_no_COVA.SNX.mat')

%% form tables: name, CRD, VEL
set_ref_sites = {'GRAS', 'GRAZ', 'LROC', 'MEDI', 'WTZR', 'ZIMM'};

% for IGS
Site = IGS_SNX.SOLUTION.ESTIMATE.StationData.Site;
SolN = IGS_SNX.SOLUTION.ESTIMATE.StationData.SolN;
CRD =  IGS_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL =  IGS_SNX.SOLUTION.ESTIMATE.Data.VEL;
NameSolN = cellstr([ cell2mat(Site), cell2mat(SolN) ]);
IGS = table(Site, CRD, VEL, 'RowNames',NameSolN);
IGS = sortrows(IGS,1);
IGS_ref = IGS(ismember(IGS.Site, set_ref_sites),:);

%% for IGB
Site = IGB_SNX.SOLUTION.ESTIMATE.StationData.Site;
SolN = IGB_SNX.SOLUTION.ESTIMATE.StationData.SolN;
CRD =  IGB_SNX.SOLUTION.ESTIMATE.Data.CRD;
VEL =  IGB_SNX.SOLUTION.ESTIMATE.Data.VEL;
NameSolN = cellstr([ cell2mat(Site), cell2mat(SolN) ]);
IGB = table(Site, CRD, VEL, 'RowNames',NameSolN);
IGB = sortrows(IGB,1);
IGB_ref = IGB(ismember(IGB.Site, set_ref_sites),:);

%% for EPN
Site = EPN_SNX.SOLUTION.ESTIMATE.StationData.Site;
SolN = EPN_SNX.SOLUTION.ESTIMATE.StationData.SolN;
CRD =  EPN_SNX.SOLUTION.ESTIMATE.Data.CRD; 
VEL =  EPN_SNX.SOLUTION.ESTIMATE.Data.VEL; 
NameSolN = cellstr([ cell2mat(Site), cell2mat(SolN) ]);
EPN = table(Site, CRD, VEL, 'RowNames',NameSolN);
EPN = sortrows(EPN,1);
EPN_ref = EPN(ismember(EPN.Site, set_ref_sites),:);

%% prepare short (average values) tables

IGB_rs = table_average(IGB_ref);
IGS_rs = table_average(IGS_ref);
EPN_rs = table_average(EPN_ref);

%% get CRD diff
% get diff in xyz and transform to -> enu

iiComm = [1:3,5:6];

[dbs_e, dbs_n, dbs_u] = XYZ2ENU((IGB_rs.CRD           + IGS_rs.CRD)          /2 ,IGB_rs.CRD           - IGS_rs.CRD) ;
[dse_e, dse_n, dse_u] = XYZ2ENU((IGS_rs.CRD(iiComm,:) + EPN_rs.CRD)          /2 ,IGS_rs.CRD(iiComm,:) - EPN_rs.CRD) ;
[deb_e, deb_n, deb_u] = XYZ2ENU((EPN_rs.CRD           + IGB_rs.CRD(iiComm,:))/2 ,EPN_rs.CRD           - IGB_rs.CRD(iiComm,:));

%%
table(IGB_rs.Site, [dbs_e, dbs_n, dbs_u]*1000)
bs = mean([dbs_e, dbs_n, dbs_u]*1000)
rms([dbs_e, dbs_n, dbs_u]*1000)

table(IGB_rs.Site(iiComm), [dse_e, dse_n, dse_u]*1000)
se = mean([dse_e, dse_n, dse_u]*1000)
rms([dse_e, dse_n, dse_u]*1000)

table(IGB_rs.Site(iiComm), [deb_e, deb_n, deb_u]*1000)
eb = mean([deb_e, deb_n, deb_u]*1000)
rms([deb_e, deb_n, deb_u]*1000)

dv_bs = mean((IGB_rs.VELenu-IGS_rs.VELenu)*1000)
dv_se = mean((IGS_rs.VELenu(iiComm,:)-EPN_rs.VELenu)*1000)
dv_eb = mean((EPN_rs.VELenu-IGB_rs.VELenu(iiComm,:))*1000)

%%
clc
fprintf('               dE       dN        du\n')
fprintf('IGB - IGS: %8.2f %8.2f %8.2f\n', bs)
fprintf('IGS - EPN: %8.2f %8.2f %8.2f\n', se)
fprintf('EPN - IGB: %8.2f %8.2f %8.2f\n', eb)
fprintf('               dVE      dVN      duV\n')
fprintf('IGB - IGS: %8.2f %8.2f %8.2f\n', dv_bs)
fprintf('IGS - EPN: %8.2f %8.2f %8.2f\n', dv_se)
fprintf('EPN - IGB: %8.2f %8.2f %8.2f\n', dv_eb)


%% get VEL in ENU

[IGB_Ve, IGB_Vn, IGB_Vu] = XYZ2ENU(IGB_rs.CRD, IGB_rs.VEL) ;
VELenu = [IGB_Ve, IGB_Vn, IGB_Vu];
IGB_rs = [IGB_rs, table(VELenu)];

[IGS_Ve, IGS_Vn, IGS_Vu] = XYZ2ENU(IGS_rs.CRD, IGS_rs.VEL) ;
VELenu = [IGS_Ve, IGS_Vn, IGS_Vu];
IGS_rs = [IGS_rs, table(VELenu)];

[EPN_Ve, EPN_Vn, EPN_Vu] = XYZ2ENU(EPN_rs.CRD, EPN_rs.VEL) ;
VELenu = [EPN_Ve, EPN_Vn, EPN_Vu];
EPN_rs = [EPN_rs, table(VELenu)];


%% compare velocities
close all
clc

clr = lines(6);
ii = [1:3,5,6]';
close all
fig1 = figure(1);
subplot(3,3,1)
hold on; grid on; axis equal; axis square
de = -(IGB_rs.VELenu(ii,1) - EPN_rs.VELenu(:,1));
xlim([18 22])
ylim([18 22])
plot([18 22], [18 22])
for i = 1:5
    plot(IGB_rs.VELenu(ii(i),1)*1000,  EPN_rs.VELenu(i,1)*1000 ,'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de*1000),'%5.3f'), ' mm/yr'])

%
subplot(3,3,2)
hold on; grid on; axis equal; axis square
dn = -(IGB_rs.VELenu([1:3,5,6],2) - EPN_rs.VELenu(:,2))*1000;
xlim([15 18])
ylim([15 18])
plot([15 18], [15 18])
for i = 1:5
    plot(IGB_rs.VELenu(ii(i),2)*1000, EPN_rs.VELenu(i,2)*1000,'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.3f'), ' mm/yr'])

subplot(3,3,3)
hold on; grid on; axis equal; axis square
du = -(IGB_rs.VELenu([1:3,5,6],3) - EPN_rs.VELenu(:,3))*1000;
xlim([-.5 1.5])
ylim([-.5 1.5])
plot([-.5 1.5], [-.5 1.5])
for i = 1:5
    plot(IGB_rs.VELenu(ii(i),:)*1000, EPN_rs.VELenu(i,3)*1000,'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('IGB')
ylabel('EPN')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(du),'%5.3f'), ' mm/yr'])

subplot(3,3,4)
hold on; grid on; axis equal; axis square
de = (IGB_rs.VELenu(:,1) - IGS_rs.VELenu(:,1))*1000;
xlim([18 22.5])
ylim([18 22.5])
plot([18 22.5], [18 22.5])
for i = 1:6
    plot(IGB_rs.VELenu(i,1)*1000, IGS_rs.VELenu(i,1)*1000,'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:))
end
xlabel('IGB')
ylabel('IGS')
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de),'%5.3f'), ' mm/yr'])

subplot(3,3,5)
hold on; grid on; axis equal; axis square
dn = (IGB_rs.VELenu(:,2) - IGS_rs.VELenu(:,2))*1000;
xlim([15 18])
ylim([15 18])
plot([15 18], [15 18])
for i = 1:6
    plot(IGB_rs.VELenu(i,2)*1000, IGS_rs.VELenu(i,2)*1000,'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:))
end
xlabel('IGB')
ylabel('IGS')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.3f'), ' mm/yr'])

s6 = subplot(3,3,6);
hold on; grid on; axis equal; axis square
dn = (IGB_rs.VELenu(:,3) - IGS_rs.VELenu(:,3))*1000;
xlim([-1.5 1.5])
ylim([-1.5 1.5])
plot([-1.5 1.5], [-1.5 1.5])
plt = zeros(6,1);
for i = 1:6
    plt(i) = plot(IGB_rs.VELenu(i,3)*1000, IGS_rs.VELenu(i,3)*1000,'o',  'MarkerFaceColor',clr(i,:), 'MarkerEdgeColor',clr(i,:));
end
xlabel('IGB')
ylabel('IGS')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(dn),'%5.3f'), ' mm/yr'])

subplot(3,3,7)
hold on; grid on; axis equal; axis square
de = (IGS_rs.VELenu([1:3,5,6],1) - EPN_rs.VELenu(:,1))*1000;
xlim([18 22.5])
ylim([18 22.5])
plot([18 22.5], [18 22.5])
for i = 1:5
    plot(EPN_rs.VELenu(i,1)*1000, IGS_rs.VELenu(ii(i),1)*1000, 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_e [mm/yr], mean_d_V_e = ', num2str(mean(de),'%5.3f'), ' mm/yr'])

subplot(3,3,8)
hold on; grid on; axis equal; axis square
dn = (IGS_rs.VELenu([1:3,5,6],2) - EPN_rs.VELenu(:,2))*1000;
xlim([15 18])
ylim([15 18])
plot([15 18], [15 18])
for i = 1:5
    plot(EPN_rs.VELenu(i,2)*1000, IGS_rs.VELenu(ii(i),2)*1000, 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_n [mm/yr], mean_d_V_n = ', num2str(mean(dn),'%5.2f'), ' mm/yr'])

s9 = subplot(3,3,9);
hold on; grid on; axis equal; axis square
du = (IGS_rs.VELenu([1:3,5,6],3) - EPN_rs.VELenu(:,3))*1000;
xlim([-.5 1.5])
ylim([-.5 1.5])
plot([-.5 1.5], [-.5 1.5])
for i = 1:5
    plot(EPN_rs.VELenu(i,3)*1000, IGS_rs.VELenu(ii(i),3)*1000, 'o',  'MarkerFaceColor',clr(ii(i),:), 'MarkerEdgeColor',clr(ii(i),:))
end
xlabel('EPN')
ylabel('IGS')
title(['V_u [mm/yr], mean_d_V_u = ', num2str(mean(du),'%5.2f'), ' mm/yr'])
leg = legend(s9, plt, 'GRAS','GRAZ','LROC','MEDI','WTZR','ZIMM','location','NorthWest');
set(leg,'FontSize',6)

%%
print(fig1,'-depsc', '-r300','ComparisonOfRefSites_vel_new.eps')


