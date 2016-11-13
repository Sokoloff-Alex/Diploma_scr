% compute strain field from deformaion

close all
clear all
clc

%% load defomation
Deformation = struct2array(load('dat/Alps_deformation_0.25x0.25_no_correlaion_3.mat'));

Vel = Deformation(:,[3,4]);
LongGrid = Deformation(:,1);
LatGrid = Deformation(:,2);
%%
Lat1  = LatGrid;
Long1 = LongGrid;
Lat2  = zeros(size(Lat1));
Long2 = zeros(size(Long1));

for i = 1:length(Lat1)
    Lat2(i)   = Lat1(i) + km2deg(Vel(i,1)/1000);
    Long2(i) = Long1(i) + km2deg(Vel(i,2)/1000, 6378*cosd(Lat1(i)) );    
end

dx = Long2 - Long1; % [deg]
dy = Lat2  - Lat1 ; % [deg]

%% Deformation tensor F
clc
clear Strain lat1 long1

p = 0;
for i = 1:10 %length(Lat1)
    fxx = (Lat2(i)  - Lat1(i))  / Lat1(i);
    fyy = (Long2(i) - Long1(i)) / Long1(1);
    fyx = dx(i) / (Lat1(i)  + dy(i));
    fxy = dy(i) / (Long1(i) + dx(i));
    F = [fxx fxy ; fyx, fyy];
    if det(F) <= 0
         disp(['Error: det(F) < 0 !!! , det(F) = ', num2str(det(F))]);
    else
        %%  Rotation matrix R
        % R = 1/2 * (F - F')
        w = 1/2 * (fyx - fxy);
        displacement = deg2km(w, cosd(Lat1(i)))*1000*1000; % disp [mm]

        %%  Strain Tensor E,  E = 1/2*(F + F')  = [exx exy; eyx, eyy 
        E =  [fxx fxy+w ; fyx-w, fyy];
        exx = E(1,1);
        eyy = E(2,2);
        exy = E(1,2); % = eyx

        %% HauptStrains (lambda1 labmbda2)   
        % [exx-L1,  exy ;
        %  eyx,  eyy-L2 ]
        % solve ax^2 + bx +c = 0
        a = 1;
        b = - exx - eyy;
        c = exx*eyy - exy^2;

        L1 = ( -b + sqrt(b^2 + 4*a*c)) / 2*a;
        L2 = ( -b - sqrt(b^2 + 4*a*c)) / 2*a;
        %  check, OK
        %  e1 = sqrt( (exx - eyy)^2 + 4*exy^2 );
        %  L1c = 1/2 * (exx + eyy + e1);
        %  L2c = 1/2 * (exx + eyy - e1);      
%%%
%         exx = 8.39;
%         exy = 4.16;
%         eyy = -1.26;
%         L1 = 9.94;
%         L2 = -2.8;
%%%
        % Omega, Eigenwert 
        Omega1  = atand( -(exx - L1) / exy )
        %  Omega2 = atand( -(exy) / (eyy - L1) ); % omega1 + omega2 = 90 deg ! 
        Omega1c = atand( (2*exy) / (exx - eyy) ) / 2   % ERROR, Omega1c ~= Omega1 

        if min(isreal([Omega1, L1, L2]))
            p = p + 1;
            Strain(p,:) = [Omega1, L1, L2];
            lat1(p,1) = Lat1(i);
            long1(p,1) = Long1(i);
        end
	end
end
% Strain;
disp('done')

%%  
plotStrain(Strain, long1, lat1, 100*1000*1000, Z, refvec);








