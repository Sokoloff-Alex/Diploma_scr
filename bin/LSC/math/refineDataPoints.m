function [x, y] = refineDataPoints(x,y)
% function to add more artificial points b/w 0 and 0..1 classes
%
% Input     :   x,y - initial data points
% Output    :   x,y - data points with artificial set
% 
% command: 
%   [x, y] = refineDataPoints(x,y)
%
% Alexandr sokolov, KEG
% 27.10.2016

%% approximate 1st segment as linear f(x) = a*x + b
b = y(1) - y(2);
a = (y(1) - y(2)) / (x(1) - x(2));

% make artificial points
step = (x(2) - x(1))/4;

x1 = (x(1)+step:step:x(2)-step)';

y1 = a.*x1 + b;

x = [x(1); x1; x(2:end)];
y = [y(1); y1; y(2:end)];

end
