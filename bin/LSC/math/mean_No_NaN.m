function [x, y] = mean_No_NaN(Corr_Mat)
% trim columns with NaN and compute mean columnwise
% since NaNs produce artificial data points
% 
% Corr_Mat  - Matrix with correlations 
% y         - vector mean values columnwise
% x         - columns used, started from 0

% Example
%   [x, y] = mean_No_NaN(Corr_Mat)
%
% Alexandr Sokolov, KEG
% 28.10.2016

y1 = nanmean(Corr_Mat);
x1 = 0:1:size(Corr_Mat,2)-1;

x = x1(~isnan(y1));
y = y1(~isnan(y1));

% p = size(Corr_Mat,1);
% i = 0;
% for column = 1:size(Corr_Mat,2)
%    if ~min(isnan(Corr_Mat(:,column))) 
%       i = i + 1;
%       x(i,1) = (column-1);
%       y(i,1) = nanmean(Corr_Mat(:,column))/p;
%    end
% end

end