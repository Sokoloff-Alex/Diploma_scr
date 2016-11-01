function [V_pred] = solve_LSC_3D(lat0, long0, lat, long, Ve, Vn, Vu, R_aoi)        
% solve Least Square Collocation
%
% input : lat0, long0 - coordinates of grid point , [deg]
%         lat, long   - coordinates of observation points
%         Ve Vn Vu    - velocity components
%         R_aoi       - radius of AOI, [km]  
% output  V_pred      - interpolated / predicted velocity vector, ENU
%
%
% Alexandr Sokolov, KEG
% 19.10.2016

        
%% COMPUTE CORRELATIONS
% for 7 classes, 0 (d=0), 0..1, 1..2, 2..3, .. 5..6
% prepare correlations

Vn_res = Vn - mean(Vn);
Ve_res = Ve - mean(Ve);
Vu_res = Vu - mean(Vu);

p = length(lat);
nClasses = 7;
Corr_NN = NaN(p,nClasses);
Corr_EE = NaN(p,nClasses);
Corr_UU = NaN(p,nClasses);
Corr_NE = NaN(p,nClasses);
Corr_NU = NaN(p,nClasses);
Corr_EU = NaN(p,nClasses);

arc2 = zeros(p);

scale = R_aoi*2/(nClasses-2); % since max dist in R_aoi = diameter

for i = 1:p
    arc = distance(lat(i), long(i), lat(i:p), long(i:p))*111 ; % km
    arc2(:,i) = distance(lat(i), long(i), lat, long)*111 ; % km
    for c = 0:nClasses-1
        J = i:p;
        sel1 = J(arc <= c*scale);    % sel1
        sel2 = J(arc >  (c-1)*scale); % sel2
        J = intersect(sel1, sel2);
        n   = length(J);
        if n > 0
            for j = 1:length(J)
                Corr_NN(i+j-1,c+1) = Vn_res(i)*Vn_res(J(j));
                Corr_EE(i+j-1,c+1) = Ve_res(i)*Ve_res(J(j));
                Corr_NE(i+j-1,c+1) = Vn_res(i)*Ve_res(J(j));
                Corr_UU(i+j-1,c+1) = Vu_res(i)*Vu_res(J(j));
                Corr_NU(i+j-1,c+1) = Vn_res(i)*Vu_res(J(j));
                Corr_EU(i+j-1,c+1) = Ve_res(i)*Vu_res(J(j));
            end
        end
    end
end

% arc2./scale
% ceil(triu(arc2./scale))

% Corr_NN

cov_NN_values = nansum(Corr_NN)/p;
cov_EE_values = nansum(Corr_EE)/p;
cov_NE_values = nansum(Corr_NE)/p;
cov_UU_values = nansum(Corr_UU)/p;
cov_NU_values = nansum(Corr_NU)/p;
cov_EU_values = nansum(Corr_EU)/p;
Classes = 0:nClasses-1;
coeff_NN = fitExp([1,-0.1], Classes*scale, cov_NN_values);
coeff_EE = fitExp([1,-0.1], Classes*scale, cov_EE_values);
coeff_NE = fitExp([1,-0.1], Classes*scale, cov_NE_values);
coeff_UU = fitExp([1,-0.1], Classes*scale, cov_UU_values);
coeff_EU = fitExp([1,-0.1], Classes*scale, cov_EU_values);
coeff_NU = fitExp([1,-0.1], Classes*scale, cov_NU_values);

%% plot Covariance functions, y=a*exp(b*d)
% close all
% figure(1)
% hold on
% grid on
% title(['point :', num2str(lat0), ' ' num2str(long0)])
% pl1 = plot(Classes*scale, cov_NN_values, 'o-b');
% pl2 = plot(Classes*scale, cov_EE_values, 'o-r');
% pl3 = plot(Classes*scale, cov_NE_values, 'o-g');
% pl4 = plot(Classes*scale, cov_UU_values, 'o-b');
% pl5 = plot(Classes*scale, cov_EU_values, 'o-r');
% pl6 = plot(Classes*scale, cov_NU_values, 'o-g');
% d = (0:0.1:nClasses-1)*scale;
% plot(d, coeff_NN(1).*exp(coeff_NN(2)*d),'--b');
% plot(d, coeff_EE(1).*exp(coeff_EE(2)*d),'--r');
% plot(d, coeff_NE(1).*exp(coeff_NE(2)*d),'--g');
% plot(d, coeff_UU(1).*exp(coeff_UU(2)*d),'--b');
% plot(d, coeff_EU(1).*exp(coeff_EU(2)*d),'--r');
% plot(d, coeff_NU(1).*exp(coeff_NU(2)*d),'--g');
% legend([pl1 pl2 pl3 pl5 pl5 pl6],'C_N_N','C_E_E','C_N_E','C_U_U','C_E_U','C_N_U')
% hold off

%% Build Covariance matrices

C_obs_NN = zeros(p);
C_obs_EE = zeros(p);
C_obs_NE = zeros(p);
C_obs_UU = zeros(p);
C_obs_EU = zeros(p);
C_obs_NU = zeros(p);


C_new_NN = zeros(p,1);
C_new_EE = zeros(p,1);
C_new_NE = zeros(p,1);
C_new_UU = zeros(p,1);
C_new_EU = zeros(p,1);
C_new_NU = zeros(p,1);

for i = 1:p
   for j = 1:p
       d = distance(lat(i), long(i), lat(j), long(j)) * 111 ; % km
       C_obs_NN(i,j) = coeff_NN(1).*exp(coeff_NN(2)*d/100);
       C_obs_EE(i,j) = coeff_EE(1).*exp(coeff_EE(2)*d/100);
       C_obs_NE(i,j) = coeff_NE(1).*exp(coeff_NE(2)*d/100);
       C_obs_UU(i,j) = coeff_UU(1).*exp(coeff_UU(2)*d/100);
       C_obs_EU(i,j) = coeff_EU(1).*exp(coeff_EU(2)*d/100);
       C_obs_NU(i,j) = coeff_NU(1).*exp(coeff_NU(2)*d/100);
   end
   
   d = distance(lat0, long0, lat(i), long(i)) * 111 ; % km
   C_new_NN(i,1) = coeff_NN(1).*exp(coeff_NN(2)*d/100);
   C_new_EE(i,1) = coeff_EE(1).*exp(coeff_EE(2)*d/100);
   C_new_NE(i,1) = coeff_NE(1).*exp(coeff_NE(2)*d/100);
   C_new_UU(i,1) = coeff_UU(1).*exp(coeff_UU(2)*d/100);
   C_new_EU(i,1) = coeff_EU(1).*exp(coeff_EU(2)*d/100);
   C_new_NU(i,1) = coeff_NU(1).*exp(coeff_NU(2)*d/100);
end

C_obs = [C_obs_EE, C_obs_NE, C_obs_EU; ...
         C_obs_NE, C_obs_NN, C_obs_NU; ...
         C_obs_EU, C_obs_NU, C_obs_UU];

C_new = [C_new_EE, C_new_NE, C_new_EU; ...
         C_new_NE, C_new_NN, C_new_NU; ...
         C_new_EU, C_new_NU, C_new_UU];

%% Solve LSC
% observations  
V_obs = [Ve; Vn; Vu];


V_pred = C_new' * C_obs^-1 * V_obs;
        
        
end