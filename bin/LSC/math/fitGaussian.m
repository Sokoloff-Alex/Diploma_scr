function [bestx] = fitGaussian(apriori, tdata, ydata)
% function to fit normal (Gaussian) curve to data,
%
% find coeff. a and b in y = a*exp(b*x^2)
% intup  : x, y
% output : [a, b] 

fun = @(x)sseGaussian(x,tdata,ydata);
bestx = fminsearch(fun,apriori);

end