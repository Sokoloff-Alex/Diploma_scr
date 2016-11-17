% Srcipt to interpolate SNX std values

close all
clear all
clc

%%
SINEX = readSNX('STA/FMC_IGB_W7.SNX');
Lat = SINEX.SITE.ID.LAT;
Lon = SINEX.SITE.ID.LON;
CRD = SINEX.SOLUTION.ESTIMATE.Data.CRD;
CRD_std = SINEX.SOLUTION.ESTIMATE.Data.CRD_STD;
Sites_list = SINEX.SITE.ID.CODE;
Stations = cellstr(SINEX.SITE.ID.CODE);
DOMES    = cellstr(SINEX.SITE.ID.DOMES);
SiteDome = [SINEX.SITE.ID.CODE, repmat(char(' '),297,1), SINEX.SITE.ID.DOMES];
SiteDome_list = cellstr(SiteDome);

%%


