function [C] = empiricalCovariance(functionName,coeffs, d)
% compute value for empirical covariance function
% approximated by functionName 
%
% Input     : funtionName -  either 'exp1', 'Hirvonen' or 'normal' ('exp2')
%             coeffs      -  coeffitients [C0, a] for function
%             d           -  distance, [km];
% Output      C           -  value of covariance at point d, C(d)
%
% command:
% [C] = empiricalCovariance(functionName,coeffs, d)
% Alexandr Sokolov, KEG
% 01.11.2016

if strcmp(functionName, 'exp1')         % y = C0 * exp(a*x)

    C = coeffs(1) .* exp(coeffs(2)*d);
            
elseif strcmp(functionName, 'Hirvonen') % y = C0 / (1 + (a/x)^2) 
        
    C = coeffs(1) ./ ( 1 + ( d/coeffs(2)).^2 ); 
        
elseif max(strcmp(functionName, {'exp2','normal'})) % y = C0 * exp(a*x^2)
    
    C = coeffs(1) .* exp( coeffs(2).*d.^2 );
        
else
    disp('Error. Use one of the following functions: exp1, hirvonen, hormal');    
end

end