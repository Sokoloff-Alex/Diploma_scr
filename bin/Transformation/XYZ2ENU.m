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
      
    for i=1:len

            [lat(i),lon(i),h(i)] = ecef2geodetic(X(i),Y(i),Z(i),referenceEllipsoid('wgs84'));

            R = [-sin(lon(i))               cos(lon(i))                0
                 -cos(lon(i))*sin(lat(i))  -sin(lon(i))*sin(lat(i))   cos(lat(i))
                  cos(lon(i))*cos(lat(i))   sin(lon(i))*cos(lat(i))   sin(lat(i))];

            Venu = ( R*[Vx(i); Vy(i); Vz(i)] )' ;
            Ve(i) = Venu(1);
            Vn(i) = Venu(2);
            Vu(i) = Venu(3);      

    %         [Ve(i),Vn(i),Vu(i)] = ecef2enu(Vx(i), Vy(i), Vz(i),lat,lon,h,referenceEllipsoid('wgs84'));  
    end
    lat = rad2deg(lat);  
    lon = rad2deg(lon);  
end