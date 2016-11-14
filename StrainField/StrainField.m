% compute strain field from deformaion

close all
clear all
clc

%% load defomation
Deformation = struct2array(load('dat/Alps_deformation_0.25x0.25_no_correlaion_3.mat'));
% [Z, refvec] = etopo('MAP/etopo1_bed_c_f4/etopo1_bed_c_f4.flt', 1, [40 54], [-7 19]); % ETOPO

Vel = Deformation(:,[3,4]);
LongGrid = Deformation(:,1);
LatGrid  = Deformation(:,2);

%%
close all
s = 100;
Strain = getStrain(Deformation);
figure(1)
hold on; grid on; axis equal
quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'k');
% plotStrain(Strain.NormalStrain, Strain.Grid(:,1), Strain.Grid(:,2), 100*1000*1000);
plotStrain(Strain.ShearStrain, Strain.Grid(:,1), Strain.Grid(:,2), 50*1000*1000);



%%
Lat1  = LatGrid; % y
Long1 = LongGrid;% x
N = length(Lat1);
Lat2  = zeros(N,1);
Long2 = zeros(N,1);

for i = 1:N
    Lat2(i)  = Lat1(i)  + km2deg(Vel(i,1)/1000);
    Long2(i) = Long1(i) + km2deg(Vel(i,2)/1000, 6371*cosd(Lat1(i)) );    
end

dx = Long2 - Long1; % [deg]
dy = Lat2  - Lat1;  % [deg]

%% Deformation tensor F
clc
clear Strain lat1 long1

Strain = NaN(N,3);
lat1   = NaN(N,1);
long1  = NaN(N,1);
p = 0;
for i = 1:N
%     fyy = (Long2(i) - Long1(i)) / Long1(1);
%     fxx = (Lat2(i)  - Lat1(i))  / Lat1(i);
%     fyx = dx(i) / (Lat1(i)  + dy(i));
%     fxy = dy(i) / (Long1(i) + dx(i));
    fxx = dx(i) / Long1(i);
    fyy = dy(i) / Lat1(i);
    fxy = dx(i) / (Lat1(i)  + dy(i));
    fyx = dy(i) / (Long1(i) + dx(i));
    F = [fxx fxy ; fyx, fyy];
    if det(F) <= 0
%          disp(['Error: det(F) < 0 !!! , det(F) = ', num2str(det(F))]);
    else
        %%  Rotation matrix R
        % R = 1/2 * (F - F')
        w = (fyx - fxy)/2;
        displacement = deg2km(w, cosd(Lat1(i)))*1000*1000; % disp [mm]

        %%  Strain Tensor E,  E = 1/2*(F + F')  = [exx exy; eyx, eyy]; 
        E =  [fxx fxy+w ; fyx-w, fyy];
        exx = E(1,1);
        eyy = E(2,2);
        exy = E(1,2); % = eyx

        %% HauptStrains (lambda1 labmbda2)   
        % Principal normal strain
        % [exx-L1,  exy ;
        %  eyx,  eyy-L2 ]
        % solve ax^2 + bx +c = 0
        a = 1;
        b = - exx - eyy;
        c = exx*eyy - exy^2;

        L1 = ( -b + sqrt(b^2 - 4*a*c)) / (2*a);
        L2 = ( -b - sqrt(b^2 - 4*a*c)) / (2*a);
        %  check, OK
        %  e1 = sqrt( (exx - eyy)^2 + 4*exy^2 );
        %  L1c = (exx + eyy + e1) / 2;
        %  L2c = (exx + eyy - e1) / 2;      
%%%
%         exx =  8.39;
%         exy =  4.16;
%         eyy = -1.26;
%         L1 =  9.94;
%         L2 = -2.8;
%%%
        %% Omega, Eigenwert 
        Omega1  = atand( -(exx - L1) / exy );
        %  Omega2 = atand( -(exy) / (eyy - L1) ); % omega1 + omega2 = 90 deg ! 
        Omega1c = 1/2 * atand(2*exy/(exx - eyy)); % ERROR, Omega1c ~= Omega1 
        
        p = p + 1;
        if min(isreal([L1, L2, Omega1c]))
            Strain(p,:) = [L1, L2, Omega1c];
        else
            Strain(p,:) = [NaN, NaN, NaN];
        end
        lat1(p,1)  = Lat1(i);
        long1(p,1) = Long1(i);
    end
end
% Strain;
disp('done')

%%  
close all
s = 100;  % [mm/yr]
figure(1)
hold on; grid on
axis equal
% geoshow(Z, refvec, 'DisplayType', 'texturemap');
% demcmap(Z);
% cptcmap('Europe')
% Earth_coast(2)
xlim([-2 18])
ylim([42 50])
pl0 = quiver(LongGrid, LatGrid, Vel(:,1)*s, Vel(:,2)*s, 0, 'k');
[pl1, pl2, pl3, pl4] = plotStrain(Strain, long1, lat1, 100*1000*1000);
legend([pl0, pl1, pl3], 'Defomation field', 'Strain, mjr', 'Strain mnr');
title('velocity field / deformation model / strain field') 
xlabel('Longitude, [deg]')
ylabel('Latitude, [deg]')
hold off








