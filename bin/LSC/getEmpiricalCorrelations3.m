function [Corr_EE,Corr_NN,Corr_UU,Corr_EN,Corr_EU, Corr_NU, nClasses, scale] = getEmpiricalCorrelations3(lat, long, Ve, Vn, Vu)
% COMPUTE CORRELATIONS
% for minimum 7 classes, 0 (d=0), 0..1, 1..2, 2..3, ...
% According to book: Observations and Leasts Squares, E. Mikhail, 1976
% 
% input:    lat, long   - velocity components in [deg]
%           Ve, Vn, Vu  - velocity, preferably in ragne near to 1 to avoid singularity in LSC 
% 
% Output :  Corr_EE, Corr_NN, Corr_UU - matrices (tables) , [p x nClasses] 
%           Corr_EN, Corr_EU, Corr_NU - matrices (tables) , [p x nClasses]
%           nClasses    - number of classes, minimum 7
%           scale       - size of 1 class, [km]
%
% Command: 
%           [Corr_NN, Corr_EE, Corr_NE] = getEmpiricalCorrelations(lat, long, Vn, Ve)
%
% Alexandr Sokolov, KEG
% 07.12.2016


Ve_res = Ve - mean(Ve);
Vn_res = Vn - mean(Vn);
Vu_res = Vu - mean(Vu);


% Search for max distanse, since R_roi is not actual after adding more stations
p = length(lat);
baselines = zeros(p);
for i = 1:p
    baselines(:,i) = greatcircleArc(lat(i), long(i), lat, long) * 111 ; % km
end

D_max = max(max(baselines));
D_roi = ceil(D_max/50)*50;

% Adjust scale and number of classes
nClasses = 6; % start with 5 classes
scale = D_roi/(nClasses-2); 

if scale > 50 % shrink scale if nessessary
    scale = 50; 
    nClasses = ceil(D_roi/scale)+2;
end

baselines = triu(baselines);
pairs = (ceil(baselines/scale)); 

%% Compute correlations
nCorMax = p*p; %
Corr_EE = NaN(nCorMax,nClasses);
Corr_NN = NaN(nCorMax,nClasses);
Corr_UU = NaN(nCorMax,nClasses);
Corr_EN = NaN(nCorMax,nClasses);
Corr_EU = NaN(nCorMax,nClasses);
Corr_NU = NaN(nCorMax,nClasses);


for i = 1:p
   for j = i:p
      class = pairs(i,j);
      Corr_EE((i-1)*p+j,class+1) = Ve_res(i)*Ve_res(j);
      Corr_NN((i-1)*p+j,class+1) = Vn_res(i)*Vn_res(j);
      Corr_UU((i-1)*p+j,class+1) = Vu_res(i)*Vu_res(j);
      Corr_EN((i-1)*p+j,class+1) = Ve_res(i)*Vn_res(j);
      Corr_EU((i-1)*p+j,class+1) = Ve_res(i)*Vu_res(j);
      Corr_NU((i-1)*p+j,class+1) = Vn_res(i)*Vu_res(j);
   end
end

end

