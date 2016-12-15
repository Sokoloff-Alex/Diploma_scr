function [CRD,VEL] = ENU2XYZ2(lat,lon, h, Ve, Vn, Vu) 
% transform velocity vector from ENU to XYZ


len = length(lat);

VEL = zeros(len,3);

CRD = llh2xyz(lon, lat);

    for i=1:len
%         [CRD(i,1),CRD(i,2),CRD(i,3)] = geodetic2ecef(referenceEllipsoid('wgs84'),lat(i),lon(i),h(i));
        
        R = [-sind(lon(i))    -cosd(lon(i))*sind(lat(i))   cosd(lon(i))*cosd(lat(i))
              cosd(lon(i))    -sind(lon(i))*sind(lat(i))   sind(lon(i))*cosd(lat(i))
                0                           cosd(lat(i))                sind(lat(i))];

        VEL(i,:) = ( R*[Ve(i); Vn(i); Vu(i)] )';
    end
end