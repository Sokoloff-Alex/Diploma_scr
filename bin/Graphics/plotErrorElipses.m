function plotErrorElipses(flag, CovENUorSigmaENU, long, lat, Ve_res, Vn_res, scale, conf, color)
% function tp plot Error elipses on map
%
%
% Alexandr Sokolov, KEG

% Rescale Covariance Matrix accoring to SINEX
%  VARIANCE FACTOR                     1.800168208557666
%  [m/yr]^2  -> [mm/yr]^2              1000^2
%  Bernese scale                       10-100

if ismember(flag , {'Cov'})
    for i = 1:length(long)
        Cov = extractCovariance(CovENUorSigmaENU, i, [1 2], 'split');
        Cov = Cov * 1000^2 * 1.8 * 25;
        if det(Cov) > 0
            mu = [long(i) + Ve_res(i)*scale, lat(i) + Vn_res(i)*scale];
            % plor Error ellipse
            error_ellipse(Cov, mu, conf, 1 , color) 
        end
    end    
elseif ismember(flag , {'Sig'})
    for i = 1:length(long)
        Cov = diag(abs(CovENUorSigmaENU(i,:)));
        Cov = Cov;
        if det(Cov) > 0
            mu = [long(i) + Ve_res(i)*scale, lat(i) + Vn_res(i)*scale];
            % plor Error ellipse
            error_ellipse(Cov, mu, conf, 1 , color) 
        end
    end
end

end