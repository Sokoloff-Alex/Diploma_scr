function [bestx, rms] = fitGaussian(apriori, xdata, ydata)
% function to fit normal (Gaussian) curve to data,
%
% find coeff. a and b in y = a*exp(b*x^2)
% intup  : x, y
% output : [a, b] 

fun = @(x)sseGaussian(x,xdata,ydata);
bestx = fminsearch(fun,apriori);

% rms of fitting
rms = sseGaussian(bestx,xdata,ydata) / length(xdata);

end