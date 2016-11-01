function [V_pred] = solve_LSC_alt2(lat0, long0, lat, long, Vn, Ve, R_roi, flag_plot)        
% solve Least Square Collocation
%
% input   :     iLat, iLong - coordinates of grid point , [deg]
%               lat, long   - coordinates of observation points
%               Ve Vn       - velocity components
%               R_roi       - radius of AOI, [km]  
% output  :     V_pred      - interpolated / predicted velocity vector
%
% Example :
%   
%  V_pred = solve_LSC(iLat, iLong, lat, long, Ve, Vn, R_aoi)
%
% Alexandr Sokolov, KEG
% 19.10.2016

        
%% COMPUTE CORRELATIONS

Vn_mean = mean(Vn);
Ve_mean = mean(Ve);

Vn_res = Vn - Vn_mean;
Ve_res = Ve - Ve_mean;

%% classes, 0 (d=0), 0..1, 1..2, 2..3 ...
% prepare correlations
% clc
p = length(lat);
nClasses = 7;
Class_NN = NaN(p,nClasses);
Class_EE = zeros(p,nClasses);
Class_NE = zeros(p,nClasses);

scale = R_roi*2/(nClasses-2); % since max dist in R_aoi = diameter


for i = 1:p
    arc = distance(lat(i), long(i), lat(i:p), long(i:p)) *111 ; % km
%     arc2(i,:) = distance(lat(i), long(i), lat, long) ; % km
    for c = 0:nClasses-1
        J = i:p;
        sel1 = J(arc <= c*scale);    % sel1
        sel2 = J(arc > (c-1)*scale); % sel2
        sel = intersect(sel1, sel2);
        if ~isempty(sel)
            for j = 1:length(sel)
              Class_NN(i+j-1,c+1) =  Vn_res(i)*Vn_res(sel(j));
              Class_EE(i,c+1) = Class_EE(i,c+1) + Ve_res(i)*Ve_res(sel(j));
              Class_NE(i,c+1) = Class_NE(i,c+1) + Vn_res(i)*Ve_res(sel(j));
            end
        end
    end
end
Class_NN;
cov_NN_values = nanmean(Class_NN)
cov_EE_values = mean(Class_EE);
cov_NE_values = mean(Class_NE);

Classes = 0:nClasses-1;
coeff_NN = fitExp([1,-0.01], Classes*scale, cov_NN_values)
coeff_EE = fitExp([1,-0.01], Classes*scale, cov_EE_values);
coeff_NE = fitExp([1,-0.01], Classes*scale, cov_NE_values);

d = (0:0.1:nClasses-1)*scale;


%% plot covariance function
% close all
if strcmp(flag_plot, 'with plot')
    figure
    hold on
    grid on
    title(['point : lat ', num2str(lat0), ', long  ' num2str(long0), ' # obs.: ', num2str(length(lat)) ])
    pl1 = plot(Classes*scale, cov_NN_values, 'o-b');
    pl2 = plot(Classes*scale, cov_EE_values, 'o-r');
    pl3 = plot(Classes*scale, cov_NE_values, 'o-g');
    pl4 = plot(d, coeff_NN(1).*exp(coeff_NN(2)*d),'--b');
    pl5 = plot(d, coeff_EE(1).*exp(coeff_EE(2)*d),'--r');
    pl6 = plot(d, coeff_NE(1).*exp(coeff_NE(2)*d),'--g');
    a = 200 ;
%     pl7 = plot(d, coeff_NN(1)./(1 + (d/a).^2), 'x-b');
%     pl8 = plot(d, coeff_EE(1)./(1 + (d/a).^2), 'x-r');
%     pl9 = plot(d, coeff_NE(1)./(1 + (d/a).^2), 'x-g')
    legend([pl1 pl2 pl3],['C_N_N  a = ',num2str(coeff_NN(1)), ' b = ', num2str(coeff_NN(2)) ], ...
                         ['C_E_E  a = ',num2str(coeff_EE(1)), ' b = ', num2str(coeff_EE(2)) ], ...
                         ['C_N_E  a = ',num2str(coeff_NE(1)), ' b = ', num2str(coeff_NE(2)) ]) ;
    hold off
end

%% Build Cov matrix

C_obs_NN = zeros(p);
C_obs_EE = zeros(p);
C_obs_NE = zeros(p);

C_new_NN = zeros(p,1);
C_new_EE = zeros(p,1);
C_new_NE = zeros(p,1);

for i = 1:p
   for j = 1:p
       d = distance(lat(i), long(i), lat(j), long(j)) * 111 ; % km
       C_obs_NN(i,j) = coeff_NN(1).*exp(coeff_NN(2)*d/100);
       C_obs_EE(i,j) = coeff_EE(1).*exp(coeff_EE(2)*d/100);
       C_obs_NE(i,j) = coeff_NE(1).*exp(coeff_NE(2)*d/100);
   end
   
   d = distance(lat0, long0, lat(i), long(i)) * 111 ; % km
   C_new_NN(i,1) = coeff_NN(1).*exp(coeff_NN(2)*d/100);
   C_new_EE(i,1) = coeff_EE(1).*exp(coeff_EE(2)*d/100);
   C_new_NE(i,1) = coeff_NE(1).*exp(coeff_NE(2)*d/100);
end

C_obs = [C_obs_NN, C_obs_NE; ...
         C_obs_NE, C_obs_EE];

C_new = [C_new_NN, C_new_NE; ...
         C_new_NE, C_new_EE];

%% Solve LSC
% observations  
V_obs = [Vn; Ve];


V_pred = C_new' * C_obs^-1 * V_obs;
        
        
end