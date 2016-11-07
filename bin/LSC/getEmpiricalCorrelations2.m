function [Corr_NN, Corr_EE, Corr_NE, nClasses, scale] = getEmpiricalCorrelations2(lat, long, Vn, Ve)
% COMPUTE CORRELATIONS
% for minimum 7 classes, 0 (d=0), 0..1, 1..2, 2..3, ...
% According to book: Observations and Leasts Squares, E. Mikhail, 1976
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


% Vn_res = Vn - mean(Vn);
% Ve_res = Ve - mean(Ve);

Vn_res = Vn;
Ve_res = Ve;

% Search for max distanse, since R_roi is not actual after adding more stations
p = length(lat);
baselines = zeros(p);
for i = 1:p
    baselines(:,i) = distance(lat(i), long(i), lat, long) * 111 ; % km
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
Corr_NN = NaN(nCorMax,nClasses);
Corr_EE = NaN(nCorMax,nClasses);
Corr_NE = NaN(nCorMax,nClasses);

for i = 1:p
   for j = i:p
      class = pairs(i,j);
      Corr_NN((i-1)*p+j,class+1) = Vn_res(i)*Vn_res(j);
      Corr_EE((i-1)*p+j,class+1) = Ve_res(i)*Ve_res(j);
      Corr_NE((i-1)*p+j,class+1) = Vn_res(i)*Ve_res(j);
   end
end

% disp('done')

%% mistaken

% for i = 1:p 
%     arc = distance(lat(i), long(i), lat(i:p), long(i:p))* 111 ; % km
%     for c = 0:nClasses-1
%         J = i:p;
%         sel1 = J(arc <= c*scale);     % sel1
%         sel2 = J(arc >  (c-1)*scale); % sel2
%         J = intersect(sel1, sel2);
%         if ~isempty(J)
%             for j = 1:length(J)
%                 Corr_NN(i+j-1,c+1) = Vn_res(i)*Vn_res(J(j));
%                 Corr_EE(i+j-1,c+1) = Ve_res(i)*Ve_res(J(j));
%                 Corr_NE(i+j-1,c+1) = Vn_res(i)*Ve_res(J(j));
%             end
%         end
%     end
% end


end

