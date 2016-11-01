function [V_res_xyz, V_pl] = remove_plate_motion (CRD, VEL, Euler_pole)
% Remove plate motion that defined by Euler pole
% intut  : CRD, VEL in xyz frame
%          Euler pole (lat[deg], long[deg], omega[deg/yr])
% output : V_res_xyz -  residual velocity in xyz frame
%        : V_pl      -  velocity of plate at stations
%
% Alexandr Sokolov, KEG
% 10.10.2016


lat_pl  = Euler_pole(1); % [deg]
long_pl = Euler_pole(2); % [deg]
w_pl    = Euler_pole(3); % [deg/yr]

W_pl = deg2rad(w_pl) *[cosd(lat_pl)*cosd(long_pl)
                       cosd(lat_pl)*sind(long_pl)
                       sind(lat_pl)];
len = size(CRD,1);
                   
V_pl  = zeros(len,3);
V_res_xyz = zeros(len,3);

    for i = 1:len
        V_pl(i,:)  = cross(W_pl, CRD(i,:)');   % Plate velocity
        V_res_xyz(i,:) = VEL(i,:) - V_pl(i,:); % Residual velocity
    end
end