% Validation large table, Capri...



% TableVel_Valid = VelocityTable;
% save('TableVel_Valid.mat','TableVel_Valid')
load('TableVel_Valid.mat')
T = TableVel_Valid;
clear TableVel_Valid

Tconstr = {'BOR1', 'GRAZ', 'JOZE', 'KOSG', 'LAMA', 'METS', 'ONSA','PENC', 'ZIMM', 'WTZR', 'MATE', 'GOPE'};
Tconstr2 = {'BOR1', 'GRAZ', 'JOZE', 'KOSG',  'METS', 'ZIMM', 'WTZR', 'GOPE'};

range = 1:size(T,1);
iiSetTconstr = range(ismember(T.Site, Tconstr));
iiSetTconstr2 = range(ismember(T.Site, Tconstr2));
%% Common stations
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


%% conv CRD



%% remore residual plate motion
clc
Vdiff = [T.Ve(iiDat), T.Vn(iiDat)]/1000 - [Ve_res(iiDatA), Vn_res(iiDatA)]

Euler_pole_app = [-32.856   55.104 3.9405e-09]';
ii = iiSetTconstr2;

vel = [T.Ve(ii), T.Vn(ii)]/1000;
Omega_Est = plate_motion(Euler_pole_app,'geo', [T.Longitude(ii), T.Latitude(ii)],Vdiff, 100);

%%
[r,v] = ENU2XYZ2(T.Latitude,  T.Longitude, zeros(size(T.Longitude)), T.Ve/1000, T.Vn/1000, zeros(size(T.Ve)));
r     = llh2xyz( T.Longitude, T.Latitude);

[V_res_xyz, V_pl] = remove_plate_motion(r, v, Euler_pole_app);

[Ve_res1,Vn_res1, Vu_res1] = XYZ2ENU(r,V_res_xyz);
[Ve_pl1, Vn_pl1,  Vu_pl1]  = XYZ2ENU(r,V_pl);



%%
clc
s = 0.2;
close all
figure(1);
hold on
showETOPO(ETOPO_Alps.Etopo_Europe,ETOPO_Alps.refvec_Etopo);
Earth_coast(2);
plot(Orogen_Alp(:,1),Orogen_Alp(:,2),'-*m')
quiver(long,        lat,        Ve_res *1000*s, Vn_res *1000*s,0, 'b')
% quiver(T.Longitude(IN), T.Latitude(IN), T.Ve(IN)*s,        T.Vn(IN)*s,       0, 'm')
quiver(T.Longitude(iiEUR), T.Latitude(iiEUR), T.Ve(iiEUR)*s,        T.Vn(iiEUR)*s,         0, 'Color',[0 .5 0])
quiver(T.Longitude(iiSetTconstr), T.Latitude(iiSetTconstr), T.Ve(iiSetTconstr)*s,        T.Vn(iiSetTconstr)*s,         0, 'r')
quiver(T.Longitude(iiEUR), T.Latitude(iiEUR), Ve_res1(iiEUR)*1000*s, Vn_res1(iiEUR)*1000*s,0 ,'m')
% quiver(T.Longitude, T.Latitude, Ve_pl1 *1000*s, Vn_pl1 *1000*s,0 ,'g')

% text(T.Longitude, T.Latitude, T.Site)

text(T.Longitude(iiDat), T.Latitude(iiDat), T.Site(iiDat), 'Color', [1 0 0])
% text(T.Longitude(iiSetComon), T.Latitude(iiSetComon), T.Site(iiSetComon), 'Color', [0 0 1])
% quiver(T.Longitude(iiSetTconstr),T.Latitude(iiSetTconstr),T.Ve(iiSetTconstr)*s,T.Vn(iiSetTconstr)*s,0, 'Color',[1 0 0], 'LineWidth',3)
text(T.Longitude(iiSetTconstr), T.Latitude(iiSetTconstr), T.Site(iiSetTconstr), 'Color', [1 0 0])
% plot(T.Longitude(iiEUR), T.Latitude(iiEUR),'o')

xlim([-5 22])
ylim([ 40 53])
grid on

%%
clc
[r,v] = ENU2XYZ2(T.Latitude, T.Longitude, zeros(size(T.Longitude)), T.Ve/1000, T.Vn/1000, zeros(size(T.Ve)));
r = llh2xyz(T.Longitude, T.Latitude);
%%
clc
close all
figure(2)
hold on
Earth_coast(3)
% plot3(r(IN,1), r(IN,2), r(IN,3), '*r')
s = 1000*1000*500;
quiver3(r(iiDat,1), r(iiDat,2), r(iiDat,3), v(iiDat,1)*s, v(iiDat,2)*s, v(iiDat,3)*s,0, 'r')
text(   r(iiDat,1), r(iiDat,2), r(iiDat,3), T.Site(iiDat))
R_w2 = [5281631.2697 -1832026.4327 3060320.6845];
R_w3 = [5300739.9317 -1859409.1298 3010697.0385];
plot3([R_w(1),  -R_w(1) ]*1.2, [R_w(2), -R_w(2) ]*1.2, [R_w(3), -R_w(3) ]*1.2, 'm')
plot3([R_w2(1), -R_w2(1)]*1.2, [R_w2(2),-R_w2(2)]*1.2, [R_w2(3),-R_w2(3)]*1.2, 'k')
plot3([R_w3(1), -R_w3(1)]*1.2, [R_w3(2),-R_w3(2)]*1.2, [R_w3(3),-R_w3(3)]*1.2, 'b')


%%

Re = 6378*1000;
iWTZR  = range(ismember(T.Site, {'WTZR'}));
iWTZR2 = range2(ismember(names, {'WTZR'}));
v_mean = mean(v(iiSetTconstr2,:),1); % v(iWTZR,:) - VEL(iWTZR2,:)
r_mean = mean(r(iiSetTconstr2,:),1); % r(iWTZR,:)
e_r = r_mean ./ norm(r_mean)
e_v = v_mean ./ norm(v_mean)

e_w = cross(e_r, e_v)
R_w = e_w * Re;
disp(num2str(R_w, '%12.6f %12.6f %12.6f \n'))
W = cross(r_mean, v_mean) / norm(r_mean)^2;
w = norm(W)

[E_lat,E_lon] = ecef2llh(R_w(1),R_w(2),R_w(3));
disp(num2str([E_lat, E_lon], '%8.3f %8.3f \n'))
disp(['Euler_pole : ',num2str([E_lat, E_lon, w], '%8.3f %8.3f %7.4e \n')])
Euler_pole_m = [E_lat, E_lon, w]';




