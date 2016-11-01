function [CRD,VEL] = ENU2XYZ(lat,lon, h, Ve, Vn, Vu) 
% transform velocity vector from ENU to XYZ

lat = deg2rad(lat);
lon = deg2rad(lon);

len = length(lat);
CRD = zeros(len,3);
VEL = zeros(len,3);
    
    for i=1:len
        [CRD(i,1),CRD(i,2),CRD(i,3)] = geodetic2ecef(referenceEllipsoid('wgs84'),lat(i),lon(i),h(i));

        R = [-sin(lon(i))    -cos(lon(i))*sin(lat(i))   cos(lon(i))*cos(lat(i))
              cos(lon(i))    -sin(lon(i))*sin(lat(i))   sin(lon(i))*cos(lat(i))
                0                         cos(lat(i))               sin(lat(i))];

        VEL(i,:) = ( R*[Ve(i); Vn(i); Vu(i)] )';
    end
end