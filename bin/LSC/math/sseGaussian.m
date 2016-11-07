function sse = sseGaussian(coeffs, x, y)
% sum of squared error, beteen data and approximated function
   
    C0 = coeffs(1);
    b  = coeffs(2);

    sse = sum((y - C0.*exp(b*x.^2)).^2);

end