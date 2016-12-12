function [Cov_enuSNX, SigmaVenu, CorrVen, Angle] = SNX_cov_transformXYZ2ENU(SNX_cov,lat, long, flag)
% function to transform BIG covariance matrixes
% beetween coordinate frames
%
% Transform covariance matrix from XYZ to ENU
%
% NOTICE! only the main diagonal block are transformed
% todo: add transformation for off giagonal blocks
%
% flag  : CRD or VEL
%
% Alexandr Sokolov, KEG
% 02.12.2016


n = size(SNX_cov,1);
p = n/6;
% Dummies
Cov_enuSNX = zeros(3*p,3*p);
SigmaVenu  = zeros(p,3);
Angle      = zeros(p,1);
CorrVen    = zeros(p,1);

% range of values to be taket from big COVA 
if ismember(flag,{'CRD', 'crd', 'R', 'r', 'Crd'})
    j = 1:3;
elseif ismember(flag,{'VEL', 'vel', 'V', 'v', 'Vel'})
    j = 4:6;
else
    disp('wRNING, choose flag: "CRD" or "VEL" ');
    disp('use flag = VEL(default)')
    j = 4:6;
end


for i = 1:p
    
    % indices in Covariance matrix
    iiCv = (i-1)*6+j; 
    covVxyzSNX = SNX_cov([iiCv],[iiCv]);
    
    % fill NaN's in the upper triangle from lower
    covVxyzSNX(1,2) = covVxyzSNX(2,1);
    covVxyzSNX(1,3) = covVxyzSNX(3,1);
    covVxyzSNX(2,3) = covVxyzSNX(3,2);
    
    % transform covariance 
    [covVenu, SigmaVenu(i,:)] = covXYZ2ENU(covVxyzSNX, lat(i), long(i));
    ii = (i-1) + [1:3];
    Cov_enuSNX([ii],[ii]) = covVenu;
    
%     % get sigma (for GMT)
%     SigmaVenu(i,1:3) = sqrt( [covVenu(1,1), covVenu(2,2), covVenu(3,3)] );
    
    % get correlation Corr_en (for GMT)
    CorrVen(i,:) = covVenu(1,2) / (SigmaVenu(i,1) * SigmaVenu(i,2));
    
    % get angle of Semimajor axis (for GMT)
    Angle(i,:) = 90 + 1/2 * atand(2*covVenu(1,2) / (covVenu(1,1) - covVenu(2,2)));
    if SigmaVenu(1) < SigmaVenu(1);
        Angle(i,:) = Angle(i,:) + 90;
    end
end


end