function [bestx] = fitHirvonen(apriori, tdata,ydata)
% function to fit curve to data,
% see fitExp
%
% find coeff. K0 and a in y = K0 / (1 + (x/a)^2)
% intup  : x, y
% output : [K0, a] 

    fun = @(x)sseHirvonen(x,tdata,ydata);
    bestx = fminsearch(fun,apriori);
    bestx(2) = abs(bestx(2));

end