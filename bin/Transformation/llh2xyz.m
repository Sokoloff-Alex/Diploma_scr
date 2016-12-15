function [r] = llh2xyz(long, lat)

RE = 6371*1000;

r = zeros(length(long),3);

for k=1:length(long)
    
    r(k,:)=[cosd(lat(k))*cosd(long(k)), cosd(lat(k))*sind(long(k)), sind(lat(k))*(1-0.5*1/298.257223563)]*RE; 
    
end
    
    
end

