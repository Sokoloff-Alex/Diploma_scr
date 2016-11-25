function etopo_fig = showETOPO(ETOPO, RefVec)
% show ETOPO map 
%
%
% Alexandr Sokolov, KEG
% 25.11.2016

%% prepare coordinates

nLon = size(ETOPO,1);
nLat = size(ETOPO,2);

step = RefVec(1);
Lat0 = RefVec(2);
Lon0 = RefVec(3);

LonRange = [ Lon0, Lon0 + nLat/step];
LatRange = [ Lat0, Lat0 - nLon/step];

%% show ETOPO

etopo_fig = imagesc(LonRange, LatRange, flipud(ETOPO)); 
set(gca,'YDir','normal'); 
cptcmap('Europe'); 
 

end