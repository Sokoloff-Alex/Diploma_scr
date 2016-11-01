function sse = sseHirvonen(x_apriori,x,y)
% sum of squared error, beteen data and approximated function
% C(d) = K0/(1+ (d/a)^2)   - Hirvonen Covariance function

K0 = x_apriori(1);
a = x_apriori(2);
sse = sum((y - K0./(1 + (x/a).^2) ).^2);

end