function [Omega_Est, dOmega_k, DOP, Omega_Est_stack, dOmega_k_stack] = plate_motion(Omega_approx,flag, CRD, VEL, Iterations)
% iteratively compute correction for Euler pole of Techtonic Plate
% input:    Omega_approx - initial value of Euler pole [deg, deg, deg/yr];
%           CRD : [n x 3] 
%           VEL : [n x 3]
%           Iterations
% output:   Omega_Est  - estimated value of Euler pole
%           dOmega_k - correction of Euler pole
%           DOP      - cov matrix
%           Omega_Est_stack, dOmega_k_stack  stacked of values from each
%           iteration
% Alexandr Sokolov, KEG
% 10.10.2016

%%

Omega_Est_stack = zeros(Iterations,3);
dOmega_k_stack  = zeros(Iterations,3);
RE = 3678*1000;

dt = 1; % [yr]

if ismember(flag, {'ECEF'}) || size(CRD,2) == 3
    disp('Coordinates / Velocities in ECEF')
    [Ve, Vn, Vu, lat, long,  h ] = XYZ2ENU(CRD,       VEL);
    [Ve, Vn, Vu, lat2,long2, h2] = XYZ2ENU(CRD+VEL*dt,VEL);

    v_lat  = (lat2  - lat) /dt; % rate, [deg/yr]
    v_long = (long2 - long)/dt; % rate, [deg/yr]
  
elseif ismember(flag, {'geo','geodetic','geocentric','LLH'}) || size(CRD,2) == 2
    disp('Coordinates / Velocities in geo')
    long   = CRD(:,1);
    lat    = CRD(:,2);
    
    v_long = zeros(size(long));
    v_lat  = zeros(size(long));
    
    % convertr velocities from [m/yr] -> [ged/yr]
    for i = 1:size(VEL,1)
        v_long(i,1) = km2deg(VEL(i,1)/1000, cosd(lat(i)) * RE);
        v_lat(i,1)  = km2deg(VEL(i,2)/1000, RE);  
    end
%     [VEL, v_long, v_lat]
    
else
    disp('WARNING: undefinded coordinate system: ECEF of geo (Long, Lat)');
end
l = zeros(length(lat)*2,1);
A = zeros(length(lat)*2,3);

disp(['Lat = ',  num2str(Omega_approx(1),'% 12.6f'), ' [deg]; ', ...
      'Long = ', num2str(Omega_approx(2),'% 12.6f'), ' [deg]; ', ...
      'w = ',    num2str(Omega_approx(3)*10^6,'%12.6f'), ' [ged/10^6 yr]' ]);

dOmega_k = [1, 1, 1]'; % arbitrary, to pass througout while loop at 1st iter
%% Iterate
% for iter = 1:Iterations
iter = 0;
while (abs(dOmega_k(1)) > 10^-4) || (abs(dOmega_k(2)) > 10^-4) || (abs(dOmega_k(3)) > 10^-10)
    iter = iter + 1;
    lat_0_k  = Omega_approx(1);
    long_0_k = Omega_approx(2);
    w_0_k    = Omega_approx(3);

    %% construct Vector of Knowns l and Desing Matrix A
    for i = 1:length(lat)
        % vector of knowns l
        l([2*i-1, 2*i],1) = [ v_lat(i)  -  cosd(lat_0_k) * sind( long(i) - long_0_k ) * w_0_k * dt ; ... 
                              v_long(i) - (sind(lat_0_k) - cosd( long(i) - long_0_k ) * tand(lat(i)) * cosd(lat_0_k)) * w_0_k * dt];

        % construst desing matrix A
        f_lat_i = dt * [-w_0_k *  sind(lat_0_k) * sind(long(i) - long_0_k); ...
                        -w_0_k *  cosd(lat_0_k) * cosd(long(i) - long_0_k); ... 
                           1   *  cosd(lat_0_k) * sind(long(i) - long_0_k)]';

        f_long_i = dt* [ w_0_k * ( cosd(lat_0_k) + cosd(long(i) - long_0_k) * tand(lat(i)) * sind(lat_0_k) ); ...
                         w_0_k * ( 0             - sind(long(i) - long_0_k) * tand(lat(i)) * cosd(lat_0_k) ); ...
                           1   * ( sind(lat_0_k) - cosd(long(i) - long_0_k) * tand(lat(i)) * cosd(lat_0_k) )]';

        A([2*i-1,2*i],:) = [f_lat_i; f_long_i];
    end 

    %% LSE
%     disp('Correcrions:');
%     dOmega_k = (A'*A)^-1*A'*l;
    
    % try with weigth matrix
    err = km2deg(0.25/(1000*1000)); % 0.25 mm/yr -> 2.25*10-9 deg/yr
    P = err^2 * eye(length(l));
    dOmega_k = (A'*P*A)^-1*A'*P*l;

    d_lat_k   = dOmega_k(1);
    d_long_k  = dOmega_k(2);
    d_w_k     = dOmega_k(3);

    % disp(['d_lat_k  = ', num2str(d_lat_k) , ' [deg]'])
    % disp(['d_long_k = ', num2str(d_long_k), ' [deg]'])
    % disp(['d_w_k    = ', num2str(d_w_k)   , ' [ged/yr]'])
    % disp(' ')

    DOP = (A'*A)^-1;
    COV = (A'*P*A)^-1;
    precision = sqrt(diag(COV))';
%     disp(['precision:', num2str(precision(1)),'   ', num2str(precision(2)),'   ', num2str(precision(3)) ])
%     %% Observation Error: 
%     for i = 1:length(lat)
% 
%         v_lat_i = - w_0_k * dt * sind(lat_0_k) * sind(long(i) - long_0_k) * d_lat_k  ...
%                   - w_0_k * dt * cosd(lat_0_k) * cosd(long(i) - long_0_k) * d_long_k ...
%                   +         dt * cosd(lat_0_k) * sind(long(i) - long_0_k) * d_w_k    ...
%                   - (v_lat(i)  - cosd(lat_0_k) * sind(long(i) - long_0_k) * w_0_k * dt);
% 
%         v_long_i =  w_0_k * dt * (cosd(lat_0_k) + cosd(long(i) - long_0_k) * tand(lat(i)) * sind(lat_0_k)) * d_lat_k  ...
%                   + w_0_k * dt * (    0         - sind(long(i) - long_0_k) * tand(lat(i)) * cosd(lat_0_k)) * d_long_k ...
%                   +         dt * (sind(lat_0_k) - cosd(long(i) - long_0_k) * tand(lat(i)) * cosd(lat_0_k)) * d_w_k    ...
%                   - (v_long(i) - (sind(lat_0_k) - cosd(long(i) - long_0_k) * tand(lat(i)) * cosd(lat_0_k)) * w_0_k * dt);
% 
%         % disp(['Observation Error : v_lat_',num2str(i),' =', num2str(v_lat_i), ' [deg]; v_long_',num2str(i),' = ', num2str(v_long_i), ' [deg]'] );
%     end 
    Omega_approx = Omega_approx + dOmega_k; % update value     
    Omega_Est = Omega_approx;
%     disp(['Lat = ',  num2str(Omega_Est(1),'% 12.6f'), ' [deg]; ', ...
%           'Long = ', num2str(Omega_Est(2),'% 12.6f'), ' [deg]; ', ...
%           'w = ',    num2str(Omega_Est(3)*10^6,'%12.6f'), ' [ged/10^6 yr]' ]);
    Omega_Est_stack(iter,:) = Omega_Est';
    dOmega_k_stack(iter,:)  = dOmega_k';
end
iter
Omega_Est_stack = Omega_Est_stack(1:iter-1,:);
dOmega_k_stack  = dOmega_k_stack(1:iter-1,:);

disp(['---Euler pole --- : iterations :', num2str(iter)])
disp(['Lat = ', num2str(Omega_Est(1),'% 12.6f'), ' [deg] :: dLat = ', num2str(dOmega_k(1)), ' [deg]'])
disp(['Lon = ', num2str(Omega_Est(2),'% 12.6f'), ' [deg] :: dLon = ', num2str(dOmega_k(2)), ' [deg]'])
disp(['w   = ', num2str(Omega_Est(3),'% 12e'),' [ged/yr] :: dw = ', num2str(dOmega_k(3)), ' [deg/yr]'])


end

                
                
                