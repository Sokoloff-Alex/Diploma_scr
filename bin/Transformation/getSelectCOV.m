function [CovSet] = getSelectCOV(CovVenuSNX, iSet)
% function to select Covariances only for selected set of sites iSet
% as well as the variance b/w sites 

% Alexandr Sokolov, KEG
% 02.12.2016

p = length(iSet);

% Dummy
CovSet = zeros(3*p,3*p);


for i = 1:p
    i1 = (i-1) + (1:3);
    i2 = (iSet(i)-1) + (1:3);
    for  j = 1:p
        j1 = (i-1) + (1:3);
        j2 = (iSet(j)-1) + (1:3);
        CovSet(i1,j1) = CovVenuSNX(i2,j2);
    end
end


end