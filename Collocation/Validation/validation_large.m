% Validation large table, Capri...

clc
clear all
close all

%%
% TableVel_Valid = VelocityTable;
% save('TableVel_Valid.mat','TableVel_Valid')
load('TableVel_Valid.mat')
T = TableVel_Valid;
clear TableVel_Valid

%% load EPN solution
clc
EPN_SNX = readSNX('EPN_A_IGb08_no_COVA.SNX', 'no COVA');
Omega_Eur = [55.9533, -97.4134,   2.6364e-07 ]';
%% remove plare motion from SNX solution

[V_res_snx, V_pl_snx] = remove_plate_motion(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD, EPN_SNX.SOLUTION.ESTIMATE.Data.VEL, Omega_Eur);
[Ve_res_s, Vn_res_s, Vu_res_s, lat_s,lon_s] = XYZ2ENU(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD,V_res_snx);
[Ve_pl_s,  Vn_pl_s,  Vu_pl_s ]              = XYZ2ENU(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD,V_pl_snx);

names_s = EPN_SNX.SOLUTION.ESTIMATE.StationData(:,2);
range_s = 1:length(names_s);
iiSNXcom = range_s(ismember(names_s, names));
iiALPcom = range_s(ismember(names, names_s));


%%  compute difference


v_s_mean =  zeros(length(iiALPcom),5);
range_c = 1:length(names_s);
for i = 1:length(iiALPcom)
   iiCom = range_c(ismember( names_s, names(iiALPcom(i)) ) );
   v_s_mean(i,:) = [long(iiALPcom(i)), lat(iiALPcom(i)), mean(Ve_res_s(iiCom)), mean(Vn_res_s(iiCom)), mean(Vu_res_s(iiCom))]; 
end

v_err = [Ve_res(iiALPcom), Vn_res(iiALPcom) , Vu_res(iiALPcom)] - v_s_mean(:,[3:5]);
dv = [long(iiALPcom), lat(iiALPcom), v_err ];

e_v_h = sqrt(sum(dv(:,3).^2 + dv(:,4).^2, 2))*1000; % [mm/yr]
v_err = v_err * 1000;
c_v_err = cov(v_err);

%% run collocation
[LongGrid,   LatGrid,   V_def,   rmsFit ,  V_SigPred  ] = run_Collocation(long(iiALPcom), lat(iiALPcom), [Ve_res(iiALPcom), Vn_res(iiALPcom), Vu_res(iiALPcom)], zeros(51*3,51*3), [-4 17], [42 50 ], 0.5, 200, 10);

%%
[LongGrid_s, LatGrid_s, V_def_s, rmsFit_s, V_SigPred_s] = run_Collocation(v_s_mean(:,1),  v_s_mean(:,2), v_s_mean(:,3:5), zeros(51*3,51*3), [-4 17], [42 50 ], 0.5, 200, 10);

%%

close all
fig2 = figure(2);
subplot(2,2,1)
hold on
axis equal
axis vis3d
grid on
plot3(v_err(:,1),v_err(:,2),v_err(:,3),'.r','MarkerSize',16)
plot3(v_err(:,1),v_err(:,2),-2*ones(51,1),'.b')
plot3(v_err(:,1),2*ones(51,1),v_err(:,3),'.b')
plot3(2*ones(51,1),v_err(:,2),v_err(:,3),'.b')
% text(v_err(:,1),v_err(:,2),v_err(:,3), names(iiALPcom))
plot3(mean(v_err(:,1)),mean(v_err(:,2)),mean(v_err(:,3)),'xk')
% legend(['ALPNET vs EPN, sigmaVe = ', num2str(sqrt(c_v_err(1,1))), ...
%        ' ; sigmaVn = ', num2str(sqrt(c_v_err(2,2)), '%3.2f') , ...
%        ' ; sigmaVu = ', num2str(sqrt(c_v_err(3,3)), '%3.2f') , ' [mm/yr]']);
% title(['ALPNET vs EPN'])
% error_ellipse(c_v_err, mean(v_err), 0.95, 'r') % 2 sigma, 95 % confidence
% error_ellipse(c_v_err, mean(v_err), 0.68, 'b') % 2 sigma, 95 % confidence
xlabel('East, [mm/yr]')
ylabel('North, [mm/yr]')
zlabel('Up, [mm/yr]')
xlim([-1.5 1.5])
ylim([-1.5 1.5])
zlim([-2 2.3])
view(-50,25)
hold off 

subplot(2,2,2)
hold on
axis equal
grid on
plot(v_err(:,1),v_err(:,2), '.b')
plot(mean(v_err(:,1)),mean(v_err(:,2)),'*r')
error_ellipse(cov(v_err(:,[1,2])), mean(v_err(:,[1,2])), 0.95, '--r') % 2 sigma, 95 % confidence
error_ellipse(cov(v_err(:,[1,2])), mean(v_err(:,[1,2])), 0.68, '--k') % 2 sigma, 95 % confidence
xlabel('East, [mm/yr]')
ylabel('North, [mm/yr]')
xlim([-1.1 0.5])

subplot(2,2,3)
hold on
axis equal
grid on
plot(v_err(:,1),v_err(:,3), '.b')
plot(mean(v_err(:,1)),mean(v_err(:,3)),'*r')
error_ellipse(cov(v_err(:,[1,3])), mean(v_err(:,[1,3])), 0.95, '--r') % 2 sigma, 95 % confidence
error_ellipse(cov(v_err(:,[1,3])), mean(v_err(:,[1,3])), 0.68, '--k') % 2 sigma, 95 % confidence
xlabel('East, [mm/yr]')
ylabel('Up, [mm/yr]')

subplot(2,2,4)
hold on
axis equal
grid on
pl1 = plot(v_err(:,2),v_err(:,3), '.b')
plot(mean(v_err(:,2)),mean(v_err(:,3)),'*r')
pl2 = error_ellipse(cov(v_err(:,[2,3])), mean(v_err(:,[2,3])), 0.95, '--r') % 2 sigma, 95 % confidence
pl3 = error_ellipse(cov(v_err(:,[2,3])), mean(v_err(:,[2,3])), 0.68, '--k') % 2 sigma, 95 % confidence
xlabel('North, [mm/yr]')
ylabel('Up, [mm/yr]')
legend([pl1 pl2 pl3], 'observations','2 sigma','1 sigma')

print(fig2, 'ALP_NET_vs_EPN.eps','-depsc','-r300');


%%
clc
disp(num2str([mean(v_err)], '%8.3f %8.3f %8.3f \n' ))
disp(num2str([rms(v_err)],  '%8.3f %8.3f %8.3f \n' ))
disp(num2str([std(v_err)],  '%8.3f %8.3f %8.3f \n' ))


%% or merge solutions

[CRD_s_ave,VEL_s_ave,uniq_names_s] = merge_stations(EPN_SNX.SOLUTION.ESTIMATE.Data.CRD, EPN_SNX.SOLUTION.ESTIMATE.Data.VEL,names_s);
[V_res_snx, V_pl_snx] = remove_plate_motion(CRD_s_ave,VEL_s_ave, Omega_Eur);

names_s = uniq_names_s;
range_s = 1:length(names_s);
iiSNXcom = range_s(ismember(names_s, names));
iiALPcom = range_s(ismember(names, names_s));



%%
Tconstr  = {'BOR1', 'GRAZ', 'JOZE', 'KOSG', 'LAMA', 'METS', 'ONSA','PENC', 'ZIMM', 'WTZR', 'MATE', 'GOPE'};
Tconstr2 = {'BOR1', 'GRAZ', 'JOZE', 'KOSG', 'METS', 'ZIMM', 'WTZR', 'GOPE'};

range = 1:size(T,1);
iiSetTconstr  = range(ismember(T.Site, Tconstr));
iiSetTconstr2 = range(ismember(T.Site, Tconstr2));
% Common stations
iiSetComon = range(ismember(T.Site, names));

%% allocate stable part in 2nd velovity field comparable to 1st

X = [-4  -4  18 18 10.4  8.9  5.3   6.3 -4];
Y = [44  54  54 50 48.5 48.5 46.5  40.0 44];
IN = range(inpolygon(T.Longitude,T.Latitude,X,Y));
iiOut2 = range(ismember(T.Site, {'DRES', 'VFCH', 'BSCN', 'KARL', 'POUS', 'PLYM', 'ESCO','LLIV','CREU','CANT'}));
IN = setdiff(IN, iiOut2);

iiDat = range(ismember(T.Site, {'WTZR','ZIMM'}));
range2 = 1:length(lat);
iiDatA = range2(ismember(names,{'WTZR','ZIMM'}));

iiEUR   = range(ismember(T.EUR,  'X'));
iiCEGNR = range(ismember(T.CEGNR,'X'));
iiAMO   = range(ismember(T.AMO,  'X'));
iiARE   = range(ismember(T.ARE,  'X'));
iiSK    = range(ismember(T.SK,   'X'));
iiUPA   = range(ismember(T.UPA,  'X'));
iiYUK   = range(ismember(T.YUK,  'X'));

%% estim Eiler pole, CRUDE
Re = 6378*1000; % m
iWTZR  = range(ismember(T.Site, {'WTZR'}));
iWTZR2 = range2(ismember(names, {'WTZR'}));
% v_mean = mean(v(iiSetTconstr2,:),1); 
% r_mean = mean(r(iiSetTconstr2,:),1); 

[r,v] = ENU2XYZ2(T.Latitude, T.Longitude, zeros(size(T.Longitude)), T.Ve, T.Vn, zeros(size(T.Ve))); % m, m/yr
[r2,v2] = ENU2XYZ2(lat, long, zeros(size(lat)), Ve_res, Vn_res, Vu_res); 
r_mean = r(iWTZR,:);
v_mean = v(iWTZR,:)%  - v2(iWTZR2);
% v_mean = -v_mean;

[Ve1, Vn1, Vu1] = XYZ2ENU(r_mean,v_mean);


e_r = r_mean ./ norm(r_mean);
e_v = v_mean ./ norm(v_mean);
e_w = cross(e_r, e_v);
R_w = e_w * Re;
W = cross(r_mean, v_mean) / norm(r_mean)^2;
w = norm(W);

[E_lat,E_lon] = ecef2llh(R_w(1),R_w(2),R_w(3));
disp(['Euler_pole : ', num2str([E_lat, E_lon, w], '%8.3f %8.3f %7.4e \n')])
Euler_pole = [E_lat, E_lon, w]';

%% improve Euler pole

Euler_pole_app = [-32.856 55.104 3.9405e-09]';
ii = iiDat;
Vdiff = [T.Ve(iiDat), T.Vn(iiDat)] - [Ve_res(iiDatA), Vn_res(iiDatA)];

% vel = [T.Ve(ii), T.Vn(ii)] - [Ve_res(iWTZR2), Vn_res(iWTZR2)];
Omega_Est = plate_motion(Euler_pole_app,'geo', [T.Longitude(ii), T.Latitude(ii)],Vdiff, 100);

%% remove residual trend

[V_res1, V_pl1] = remove_plate_motion(r, v, [Euler_pole([1,2]); Euler_pole(3)*50]);

[Ve_res1, Vn_res1, Vu_res1] = XYZ2ENU(r,V_res1);
[Ve_pl1,  Vn_pl1,  Vu_pl1 ] = XYZ2ENU(r,V_pl1);


%% plot map
clc
clr = lines(8);
s = 1000*0.25;
close all
figure(1);
hold on
% showETOPO(ETOPO_Alps.Etopo_Europe,ETOPO_Alps.refvec_Etopo);
Earth_coast(2);
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'-*m')
% quiver(long,                      lat,                      Ve_res*s,             Vn_res*s,             0, 'k')
quiver(long(iiALPcom),            lat(iiALPcom),            Ve_res(iiALPcom)*s,   Vn_res(iiALPcom)*s,   0, 'b', 'LineWidth',2)
% quiver(T.Longitude(IN),           T.Latitude(IN),           T.Ve(IN)*s,           T.Vn(IN)*s,           0, 'm')

% quiver(T.Longitude(iiEUR),        T.Latitude(iiEUR),        T.Ve(iiEUR)*s,        T.Vn(iiEUR)*s,        0, 'Color',clr(1,:))
% quiver(T.Longitude(iiAMO),        T.Latitude(iiAMO),        T.Ve(iiAMO)*s,        T.Vn(iiAMO)*s,        0, 'Color',clr(2,:))
% quiver(T.Longitude(iiARE),        T.Latitude(iiARE),        T.Ve(iiARE)*s,        T.Vn(iiARE)*s,        0, 'Color',clr(3,:))
% quiver(T.Longitude(iiSK),         T.Latitude(iiSK),         T.Ve(iiSK)*s,         T.Vn(iiSK)*s,         0, 'Color',clr(4,:))
% quiver(T.Longitude(iiUPA),        T.Latitude(iiUPA),        T.Ve(iiUPA)*s,        T.Vn(iiUPA)*s,        0, 'Color',clr(5,:))
% quiver(T.Longitude(iiYUK),        T.Latitude(iiYUK),        T.Ve(iiYUK)*s,        T.Vn(iiYUK)*s,        0, 'Color',clr(6,:))
% quiver(T.Longitude(iiCEGNR),      T.Latitude(iiCEGNR),      T.Ve(iiCEGNR)*s,      T.Vn(iiCEGNR)*s,      0, 'Color',clr(7,:))
% quiver(lon_s,                     lat_s,                    Ve_res_s*s,           Vn_res_s*s,           0, 'r')
quiver(lon_s(iiSNXcom),           lat_s(iiSNXcom),          Ve_res_s(iiSNXcom)*s, Vn_res_s(iiSNXcom)*s, 0, 'color',[0 .5 0], 'LineWidth',2)
quiver(dv(:,1),                   dv(:,2),                  dv(:,3)*s,            dv(:,4)*s,            0, 'r', 'LineWidth',2)
quiver(v_s_mean(:,1),             v_s_mean(:,2),            v_s_mean(:,3)*s,      v_s_mean(:,4)*s,      0, 'g', 'LineWidth',2)
text(dv(:,1),dv(:,2)-0.2, num2str(sqrt(sum(dv(:,3).^2 + dv(:,4).^2, 2))*1000, '%3.2f' ) )
text(dv(:,1),dv(:,2), names(iiALPcom))

% quiver(T.Longitude(iiSetTconstr), T.Latitude(iiSetTconstr), T.Ve(iiSetTconstr)*s, T.Vn(iiSetTconstr)*s, 0, 'r')
% quiver(T.Longitude(iiEUR),        T.Latitude(iiEUR),        Ve_res1(iiEUR)*s,     Vn_res1(iiEUR)*s,     0 ,'m')
% quiver(T.Longitude,               T.Latitude,               Ve_pl1*s,             Vn_pl1*s,             0 ,'g')
% quiver(T.Longitude(iWTZR),        T.Latitude(iWTZR),        Ve1*s,                Vn1*s,                0, 'y')
% text(T.Longitude, T.Latitude, T.Site)
% text(T.Longitude(iiDat), T.Latitude(iiDat), T.Site(iiDat), 'Color', [1 0 0])
% text(T.Longitude(iiSetComon), T.Latitude(iiSetComon), T.Site(iiSetComon), 'Color', [0 0 1])
% quiver(T.Longitude(iiSetTconstr),T.Latitude(iiSetTconstr),T.Ve(iiSetTconstr)*s,T.Vn(iiSetTconstr)*s,0, 'Color',[1 0 0], 'LineWidth',3)
% text(T.Longitude(iiSetTconstr), T.Latitude(iiSetTconstr), T.Site(iiSetTconstr), 'Color', [1 0 0])
% plot(T.Longitude(iiEUR), T.Latitude(iiEUR),'o')
% text(long(iiALPcom), lat(iiALPcom), names(iiALPcom))
% quiver(LongGrid_s, LatGrid_s, V_def_s(:,1)*s, V_def_s(:,2)*s, 0, 'color',[0 .5 0])
% quiver(LongGrid,   LatGrid  , V_def(:,1)*s,   V_def(:,2)*s,  0,'b')
quiver(LongGrid,   LatGrid  , dvLSC(:,1)*s,   dvLSC(:,2)*s,   0,'r')
dvLSC = V_def - V_def_s;
xlim([-6 18])
ylim([ 40 53])
grid on


%% plot 3D
clc
close all
figure(2)
hold on
Earth_coast(3)
% plot3(r(IN,1), r(IN,2), r(IN,3), '*r')
s = 1000*1000*700;
quiver3(r(iiDat,1), r(iiDat,2), r(iiDat,3), v(iiDat,1)*s, v(iiDat,2)*s, v(iiDat,3)*s,0, 'r')
text(   r(iiDat,1), r(iiDat,2), r(iiDat,3), T.Site(iiDat))
plot3([R_w(1),  -R_w(1) ]*1.2, [R_w(2), -R_w(2) ]*1.2, [R_w(3), -R_w(3) ]*1.2, 'm', 'LineWidth',2)
plot3(R_w(1),R_w(2),R_w(3),'*m')




