function km = deg2km(deg,radius_km)
% funiction to convert degrees to km
% usefull to compute length in [km] of arc of given radius in [km]

rad = rad2deg(deg);

if nargin == 1
    radius_km = 6371; % Radius of Earth
end

km = radius_km * rad; 

end