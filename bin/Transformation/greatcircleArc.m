function [Arc] = greatcircleArc(lat1, lon1, lat2, lon2, varargin)

% Calculate great circle distance between points on a sphere using the
% Haversine Formula.  LAT1, LON1, LAT2, and LON2 are in radians.  RNG is a
% length and has the same units as the radius of the sphere, R.  (If R is
% 1, then RNG is effectively arc length in radians.)

if isempty(varargin)
    Re = 6378137.0; %  Equatorial radius (6,378.1370 km)
else
    Re = varargin{1}; 
end

lat1  = deg2rad(lat1);
lon1  = deg2rad(lon1);
lat2  = deg2rad(lat2);
lon2  = deg2rad(lon2);


a = sin((lat2-lat1)/2).^2 + cos(lat1) .* cos(lat2) .* sin((lon2-lon1)/2).^2;

% Ensure that a falls in the closed interval [0 1].
a(a < 0) = 0;
a(a > 1) = 1;

rng = Re * 2 * atan2(sqrt(a),sqrt(1 - a));

Arc = rad2deg(rng/Re);

end