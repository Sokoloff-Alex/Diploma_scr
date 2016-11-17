function[ArcLength, azim] = OrthodromeArcLength(lat1,  long1, lat2, long2)
% Compute the Arc on the orthodrome and bearing (initial azimuth) 
% between Point1(lat1 long1) and other Points(lat2, long2)
%
% Indut  :      lat1, long1  - coordinates of Point1 [deg, deg]
%               lat2, long2  - set of coorrdinates of Points [deg, deg]
%
% Ouput  :      ArcLength    - Lenth of the arc on the great circle, [deg]
%               azim         - Azimut from Point1 to other Points, [deg]
%
% Alexandr Sokolov

lat1  = deg2rad(lat1);
lon1  = deg2rad(long1);
lat2  = deg2rad(lat2);
lon2  = deg2rad(long2);


ArcLength = acosd(sin(lat1).*sin(lat2) + cos(lat1).*cos(lat2).*cos(lon2-lon1));

ArcLength = real(ArcLength); % sometimes have Im component for i = j, lat1 = lat2, long1 = long2;

azim = atan2d( cos(lat2) .* sin(lon2 - lon1),...
               cos(lat1) .* sin(lat2) - sin(lat1) .* cos(lat2) .* cos(lon2-lon1) );

% Azimuths are undefined at the poles, so we choose a convention: zero at
% the north pole and pi at the south pole.

azim(lat1 <= -90) = 0;
azim(lat2 >=  90) = 0;
azim(lat2 <= -90) = pi;
azim(lat1 >=  90) = pi;

end