function [bestx, rms] = fitHirvonen(apriori, xdata,ydata)
% function to fit curve to data,
% see fitExp
%
% find coeff. K0 and a in y = K0 / (1 + (x/a)^2)
% intup  : x, y
% output : [K0, a] 

    apriori(2) = 50;
    fun = @(x)sseHirvonen(x,xdata,ydata);
    bestx = fminsearch(fun,apriori);
    bestx(2) = abs(bestx(2));
    
    % rms of fitting
    rms = sseHirvonen(bestx,xdata,ydata) / length(xdata);

end