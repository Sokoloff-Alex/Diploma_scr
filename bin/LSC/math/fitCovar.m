function [coeffs] = fitCovar(functionName, x, y, coeffs_apr)
% approximate empirical covariance function
% by one of following: 'exp1', 'Hirvonen' or 'hormal'
% Input :   functionName - one of 'exp1', 'Hirvonen' or 'hormal'
%           x, y         - datapoints
%           coeffs_apr   - initial values of parameters
% Output    coeffs       - coeffitines [C0, a]
%
% command: 
% [coeffs] = fitCovar(functionName, x, y, C0_apr, a_apr)
%
% Alexandr Sokolov, KEG
% 01.11.2016

if ismember(functionName, {'exp1', 'Hirvonen','normal','exp2','Gaussian','gaussian'})
    
    if strcmp(functionName, 'exp1')     % y = C0 * exp(a*x)
           coeffs = fitExp(coeffs_apr, x, y); 
    end
    
    if strcmp(functionName, 'Hirvonen') % y = C0 / (1 + (a/x)^2) 
           coeffs = fitHirvonen(coeffs_apr, x, y); 
    end
    
    if max(strcmp(functionName, {'exp2','normal', 'Gaussian','gaussian'})) % y = C0 * exp(a*x^2) 
           coeffs = fitGaussian(coeffs_apr, x, y); 
    end    
    
else
    disp('Error. Use one of the following functions: exp1, hirvonen, hormal');
end


end