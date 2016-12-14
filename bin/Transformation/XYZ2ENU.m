function [Ve,Vn, Vu, lat,lon, h] = XYZ2ENU(CRD,VEL)
% transform velocity vector from XYZ to ENU

    X  = CRD(:,1);
    Y  = CRD(:,2);
    Z  = CRD(:,3);
    Vx = VEL(:,1);
    Vy = VEL(:,2);
    Vz = VEL(:,3);

    len = length(X);
    Ve  = zeros(len,1);
    Vn  = zeros(len,1);
    Vu  = zeros(len,1);
    lat = zeros(len,1);
    lon = zeros(len,1);
    h   = zeros(len,1);
      
    for i=1:len

%             [lat(i),lon(i),h(i)] = ecef2geodetic(X(i),Y(i),Z(i),referenceEllipsoid('wgs84'));
            [lat(i),lon(i),h(i)] = ecef2llh(X(i),Y(i),Z(i));

            R = [-sind(lon(i))                cosd(lon(i))                0
                 -cosd(lon(i))*sind(lat(i))  -sind(lon(i))*sind(lat(i))   cosd(lat(i))
                  cosd(lon(i))*cosd(lat(i))   sind(lon(i))*cosd(lat(i))   sind(lat(i))];

            Venu = ( R*[Vx(i); Vy(i); Vz(i)] )' ;
            Ve(i) = Venu(1);
            Vn(i) = Venu(2);
            Vu(i) = Venu(3);      

    %         [Ve(i),Vn(i),Vu(i)] = ecef2enu(Vx(i), Vy(i), Vz(i),lat,lon,h,referenceEllipsoid('wgs84'));  
    end
%     lat = rad2deg(lat);  
%     lon = rad2deg(lon);  
end