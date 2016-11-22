% read austria

clear all
close all
clc
%%

% save('Austriatable.mat','Austriatable')
% save('APOS.mat','APOS')


load('Austriatable')
load('APOS')

%% search for coordinated in APOS table

list_common = intersect(Austriatable.Site, APOS.SITE);

clear Latitude Longitude Height
Latitude = zeros(length(list_common),1);
Longitude = zeros(length(list_common),1);
Height = zeros(length(list_common),1);
Ve = zeros(length(list_common),1);
Vn = zeros(length(list_common),1);
Vu = zeros(length(list_common),1);


for iSite = 1:length(list_common)
    Site  = list_common(iSite);
    Latitude(iSite,1)  = dmsString2degrees(cell2mat(APOS.Latitude( strcmp(Site, APOS.SITE))));
    Longitude(iSite,1) = dmsString2degrees(cell2mat(APOS.Longitude(strcmp(Site, APOS.SITE))));
    Height(iSite,1)    = APOS.Height(    strcmp(Site, APOS.SITE));
    Ve(iSite,1)        = Austriatable.Ve(strcmp(Site, Austriatable.Site));
    Vn(iSite,1)        = Austriatable.Vn(strcmp(Site, Austriatable.Site));
    Vu(iSite,1)        = Austriatable.Vu(strcmp(Site, Austriatable.Site));
    
end

Table = table(Longitude, Latitude, Ve, Vn, Vu);
Table.Properties.RowNames = list_common;


 