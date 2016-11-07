function [Corr_NN, Corr_EE, Corr_NE, nClasses, scale] = getEmpiricalCorrelations(lat, long, Vn, Ve)
% COMPUTE CORRELATIONS
% for minimum 7 classes, 0 (d=0), 0..1, 1..2, 2..3, ...
% 
% input:    lat, long   - velocity components in [deg]
%           Vn, Ve      - velocity, preferably in ragne near to 1 to avoid singularity in LSC 
% 
% Output :  Corr_NN, Corr_EE, Corr_NE - matrices (tables) , [p x nClasses]
%           nClasses    - number of classes, minimum 7
%           scale       - size of 1 class, [km]
%
% Command: 
%           [Corr_NN, Corr_EE, Corr_NE] = getEmpiricalCorrelations(lat, long, Vn, Ve)
%
% Alexandr Sokolov, KEG
% 27.10.2016

%% Search for max distanse, since R_roi is not actual after adding more stations
p = length(lat);
arc2 = zeros(p);
for i = 1:p
    arc2(:,i) = distance(lat(i), long(i), lat, long) * 111 ; % km
end
D_max = max(max(arc2));
D_roi = ceil(D_max/100)*100;

%% Adjust scale and number of classes
nClasses = 4; % start with 5 classes
scale = D_roi/(nClasses-2); 

if scale > 50 % shrink scale if nessessary
    scale = 50; 
    nClasses = ceil(D_roi/scale)+2;
end

%% Compute correlations
Corr_NN = NaN(p,nClasses);
Corr_EE = NaN(p,nClasses);
Corr_NE = NaN(p,nClasses);

Vn_res = Vn - mean(Vn);
Ve_res = Ve - mean(Ve);

% Vn_res = Vn;
% % Ve_res = Ve;

for i = 1:p
    arc = distance(lat(i), long(i), lat(i:p), long(i:p))* 111 ; % km
    for c = 0:nClasses-1
        J = i:p;
        sel1 = J(arc <= c*scale);     % sel1
        sel2 = J(arc >  (c-1)*scale); % sel2
        J = intersect(sel1, sel2);
        if ~isempty(J)
            for j = 1:length(J)
                Corr_NN(i+j-1,c+1) = Vn_res(i)*Vn_res(J(j));
                Corr_EE(i+j-1,c+1) = Ve_res(i)*Ve_res(J(j));
                Corr_NE(i+j-1,c+1) = Vn_res(i)*Ve_res(J(j));
            end
        end
    end
end

