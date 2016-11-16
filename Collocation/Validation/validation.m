% script to validate SCL result with other velocity fields

%% make structure from France Supl.Data


close all
clear all
clc

%%
% data = Francevelocityfieldsupldata;
% 
% FieldNames = {'Sites','Long','Lat','Ve','Vn','Vu', ...     
%               'SigmaVe','SigmaVn','SigmaVu', ...
%               'Period','Days','years', ...
%               'kappa_E','kappa_N','kappa_U',};
% 
% FranceVELStruct = cell2struct(data(1:end,:)', FieldNames);
% 
% save('FranceVELStruct.mat','FranceVELStruct')

FranceVELStruct = load('FranceVELStruct.mat');


%%

figure(1)
hold on; grid on
plot(str2mat(FranceVELStruct.Long)), FranceVELStruct.Lat, '*')
