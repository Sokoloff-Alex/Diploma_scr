function [Degrees] = dmsString2degrees(dmsString)
% funiction for convert String DMS into float Degrees 
% 
% command:
%   [Degrees] = dmsString2degrees(dmsString)
%
% Alexandr Sokolov, KEG
% 22.11.2016

%% Parse String DMS
DMS = strsplit(dmsString);

D = sscanf(cell2mat(DMS(1)), '%d');
M = sscanf(cell2mat(DMS(2)), '%d');
S = sscanf(cell2mat(DMS(3)), '%f');

Degrees = dms2degrees([D M S]);

end