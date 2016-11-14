function [deg] = km2deg(km, radius_km)
% override funtion km2ged
% usefull to compute arc on paralleles of givel length (dependent on latitude)
%
%   Input :  km         : distande in [km]
%            radius_km  : radius of parallel in [km]         
%   Output:  deg        : arc in [deg] 
%
% Command:
%   deg = KM2DEG(km)
%
% Alexandr Sokolov, KEG
% 14.11.2016

if nargin == 1
    radius_km = 6371; % [km] Radius of Earth
end

rad = km/radius_km;

deg = rad2deg(rad);

end
