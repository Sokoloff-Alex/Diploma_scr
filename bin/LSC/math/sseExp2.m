function sse = sseExp2(coeffs,d,C)
% sum of squared error, beteen data and approximated function
   
    C0 = coeffs(1);
    b = coeffs(2);

    sse = sum((C - C0.*exp(b*d.^2)).^2);

end