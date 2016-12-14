function [r] = llh2xyz(long, lat)

RE = 6371*1000;


r = zeros(length(long),3);

long = deg2rad(long);
lat  = deg2rad(lat);
    
for k=1:length(long)
    r(k,:)=[cos(lat(k))*cos(long(k)) cos(lat(k))*sin(long(k)) sin(lat(k))*(1-0.5*1/298.257223563)]*RE; 
end
    
    
end

