function [bestx] = fitExp2(apriori, tdata, ydata)
% function to fit normal (Gaussian) curve to data,
%
% find coeff. a and b in y = a*exp(b*x^2)
% intup  : x, y
% output : [a, b] 

fun = @(x)sseExp2(x,tdata,ydata);
bestx = fminsearch(fun,apriori);

end