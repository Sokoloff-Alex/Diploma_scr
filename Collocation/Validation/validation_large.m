% Validation large table, Capri...



% TableVel_Valid = VelocityTable;
% save('TableVel_Valid.mat','TableVel_Valid')
load('TableVel_Valid.mat')

T = TableVel_Valid;


%%
s = 0.2;
close all
figure(1)
hold on
showETOPO(ETOPO_Alps.Etopo_Europe,ETOPO_Alps.refvec_Etopo)
Earth_coast(2)
quiver(T.Longitude, T.Latitude, T.Ve*s,        T.Vn*s,       0)
quiver(long,        lat,        Ve_res*1000*s, Vn_res*1000*s,0)
xlim([-10 30])
ylim([ 40 55])
grid on

