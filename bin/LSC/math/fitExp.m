function [bestx, rms] = fitExp(apriori, xdata, ydata)
% function to fit exponential curve to data,
% https://www.mathworks.com/help/matlab/math/example-curve-fitting-via-optimization.html?searchHighlight=sseval
%
% find coeff. a and b in y = a*exp(b*x)
% intup  : x, y
% output : [a, b] 

fun = @(x)sseExp1(x,xdata,ydata);
bestx = fminsearch(fun,apriori);

% rms of fitting
rms = sseExp1(bestx,xdata,ydata)/length(xdata);

end