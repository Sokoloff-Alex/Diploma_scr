% test_ETRF

Omega_app = Omega_NNR_NUVEL_1A;  

CRD_e = cell2mat( testETRF(:,2:4));
VEL_e = cell2mat( testETRF(:,5:7));

[Ve,Vn, Vu, lat3, long3, h3] = XYZ2ENU(CRD_e,VEL_e);
[Ve,Vn, Vu, lat4, long4, h4] = XYZ2ENU(CRD_e+VEL_e*dt,VEL_e);

d_lat  = lat4  - lat3;
d_long = long4 - long3;


Iterations = 200;

Omega_app_Eur_stack = zeros(Iterations+1,3);
dOmega_k_stack = zeros(Iterations+1,3);

points = points_ref_2;

%   56.9187
%  -95.0432
%    2.7004e-07
% refined
%    56.3768
%   -96.2664
%     2.6623e-07
% 
Omega_app_Eur_stack(1,:) = Omega_app;

for iter = 2:Iterations+1

[dOmega_k, DOP, Omega_Est_e] = plate_motion(Omega_app, lat3, long3, d_lat, d_long, dt);

Omega_app = Omega_Est_e;
Omega_app_Eur_stack(iter,:) = Omega_app';
dOmega_k_stack(iter,:) = dOmega_k';

end

%%

% close all
figure(1)
subplot(2,3,[1,2,4,5])
hold on
grid on
plot(Omega_app_Eur_stack(:,2), Omega_app_Eur_stack(:,1), '.r')
plot(Omega_Est(2),Omega_Est(1),'om')
plot(Omega_NNR_NUVEL_1A(2),Omega_NNR_NUVEL_1A(1),'*b')
xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
subplot(2,3,[3,6])
hold on
grid on
% plot(0,Omega_app_Eur0(3),'*m')
plot(0,Omega_NNR_NUVEL_1A(3),'*b')
plot(1:iter,Omega_app_Eur_stack(:,3),'.b')

%%

try     
    close (plResVel) 
end
plResVel = figure(6);
hold on 
grid on
Earth_coast(2)
xlim([-6 18])
ylim([41 53])
plot(long(range_flag), lat(range_flag), '.b')
plot(long(points), lat(points), '*r','lineWidth',3)
% text(long(points), lat(points), names(points),'fontsize',16)
% quiver(long, lat, d_long*s, d_lat*s, 0)
quiver(long(range_flag), lat(range_flag), v_long_res(range_flag)*s, v_lat_res(range_flag)*s, 0,'b')
quiver(long(points), lat(points), v_long_res(points)*s, v_lat_res(points)*s, 0, 'r')
% quiver(long, lat, d_long_pl*s, d_lat_pl*s,0)
legend('Earth Coast','all stations','selected stations',  'Residual velocity')

hold off

