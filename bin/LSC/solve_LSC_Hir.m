function [V_pred] = solve_LSC_Hir(lat0, long0, lat, long, Vn, Ve, varargin)
% solve Least Square Collocation
%
% input   :     lat0, long0 - coordinates of grid point , [deg]
%               lat, long   - coordinates of observation points
%               Ve Vn       - velocity components 
%               varargin    - 'plot' to plot covariance functions
%                           - '-v' verbose mode           
% output  :     V_pred      - interpolated / predicted velocity vector
%
% Example :
%   
% V_pred = solve_LSC(iLat, iLong, lat, long, Ve, Vn, R_aoi)
%
% Alexandr Sokolov, KEG
% 25.10.2016

flags = varargin(:);
p = length(lat);

%% COMPUTE CORRELATIONS
[Corr_NN, Corr_EE, Corr_NE, nClasses, scale] = getEmpiricalCorrelations(lat, long, Vn, Ve);

%% get mean and trim columns with all NaN
[x_NN, y_NN] = mean_No_NaN(Corr_NN);
[x_EE, y_EE] = mean_No_NaN(Corr_EE);
[x_NE, y_NE] = mean_No_NaN(Corr_NE);

x_NN = x_NN*scale;
x_EE = x_EE*scale;
x_NE = x_NE*scale;

%% refine data points in 1st linear segment
if ismember('refine', flags)
    [x_NN, y_NN] = refineDataPoints(x_NN,y_NN);
    [x_EE, y_EE] = refineDataPoints(x_EE,y_EE);
    [x_NE, y_NE] = refineDataPoints(x_NE,y_NE);
end

if y_NE(1) < 0
    b_NE =  -abs(max(y_NE));
else
    b_NE =   abs(min(y_NE));
end

%% add values to tail for better fitting
if ismember('tail', flags)
    n = 1;
    x_NN = [x_NN; [x_NN(end)+scale:scale:x_NN(end)+scale*n]'];
    x_EE = [x_EE; [x_EE(end)+scale:scale:x_EE(end)+scale*n]'];
    x_NE = [x_NE; [x_NE(end)+scale:scale:x_NE(end)+scale*n]'];
    
    y_NN = [y_NN; min(y_NN)*ones(n,1)];
    y_EE = [y_EE; min(y_EE)*ones(n,1)];
    y_NE = [y_NE;     -b_NE*ones(n,1)];
end

if ismember('tail 0', flags)
    n = 5;
    x_NN = [x_NN; [x_NN(end)+scale:scale:x_NN(end)+scale*n]'];
    x_EE = [x_EE; [x_EE(end)+scale:scale:x_EE(end)+scale*n]'];
    x_NE = [x_NE; [x_NE(end)+scale:scale:x_NE(end)+scale*n]'];
    
    y_NN = [y_NN; zeros(n,1)];
    y_EE = [y_EE; zeros(n,1)];
    y_NE = [y_NE; zeros(n,1)];
end

%% fit curves

[K0_NN,  a_NN] = fitHirvonen([y_NN(1),100], x_NN, y_NN);
[K0_EE,  a_EE] = fitHirvonen([y_EE(1),100], x_EE, y_EE);
[K0_NE,  a_NE] = fitHirvonen([y_NE(1),100], x_NE, y_NE);

%% set minimum value for coeff a
if ismember('set a', flags)
    a_min = 50;
    if a_NN < a_min; a_NN = a_min; end
    if a_EE < a_min; a_EE = a_min; end
    if a_NE < a_min; a_NE = a_min; end
end

%% Output statistics
% verbose mode
if ismember('-v', flags)
    disp(['point ::', ... 
          ' lat = ',num2str(lat0, '%5.2f'), ...
          '; long = ',  num2str(long0,'%5.2f'), ...
          ' # obs.: ', num2str(p), ...
          ' :: Cnn: a = ', num2str(a_NN, '%#7.4f'), ...
          ' :: Cee: a = ', num2str(a_EE, '%#7.4f'), ...
          ' :: Cne: a = ', num2str(a_NE, '%#7.4f') ]);
end

%% plot Covariance functions, y=K0/(1 + (d/a)^2)
d = (0:(max(x_NN)/100):max(x_NN));
if ismember('plot', flags)
    figure
    hold on
    grid on
    title(['point : lat ', num2str(lat0), ', long  ' num2str(long0), ...
        ' # obs.: ', num2str(length(lat)), ' ; step = ', num2str(scale), ...
        ' km', ';  # Classes ', num2str(nClasses)])
    pl1 = plot(x_NN, y_NN, 'o--b');
    pl2 = plot(x_EE, y_EE, 'o--r');
    pl3 = plot(x_NE, y_NE, 'o--g');  
    pl4 = plot(d,K0_NN./(1+(d/a_NN).^2), '-b', 'LineWidth',2);
    pl5 = plot(d,K0_EE./(1+(d/a_EE).^2), '-r', 'LineWidth',2);
    pl6 = plot(d,K0_NE./(1+(d/a_NE).^2), '-g', 'LineWidth',2);
    pl7 = plot(d, K0_NN(1).*exp(-0.005*d ), '.-k');
    pl8 = plot(d, K0_NN(1).*exp(-0.05*d), '--k');
    legend([pl4 pl5 pl6 pl7 pl8],['C_N_N  K0 = ', num2str(K0_NN,'%.2g'), ' a = ', num2str(a_NN,'%4.0f') ], ...
                         ['C_E_E  K0 = ', num2str(K0_EE,'%.2g'), ' a = ', num2str(a_EE,'%4.0f') ], ...
                         ['C_N_E  K0 = ', num2str(K0_NE,'%.2g'), ' a = ', num2str(a_NE,'%4.0f') ], ...
                         ['stable : ', num2str(K0_NN,'%.2g'), ' *exp(',num2str(-0.005,'%.2g'), '*d)'], ...
                         ['deform : ', num2str(K0_NN,'%.2g'), ' *exp(',num2str(-0.05, '%.2g'), '*d)']);
     xlabel('Classes, (distance [km])')
     ylabel('Covariance')
     hold off
end

%% histogram of distances
if ismember('hist', flags)
    figure(2)
    hold on; grid on
    hist(dist,nClasses)
    title(['Histogram of baselines, # of sites: ', num2str(p)])
    xlabel('Classes , [km]')
    ylabel(' # of pairs')
    hold off
end

%% Build Covariance matrices
C_obs_NN = zeros(p,p);
C_obs_EE = zeros(p,p);
C_obs_NE = zeros(p,p);
C_new_NN = zeros(p,1);
C_new_EE = zeros(p,1);
C_new_NE = zeros(p,1);

for i = 1:p
   for j = 1:p
       d = distance(lat(i), long(i), lat(j), long(j)) * 111 ; % km
       C_obs_NN(i,j) = K0_NN./(1 + (d/a_NN).^2);
       C_obs_EE(i,j) = K0_EE./(1 + (d/a_EE).^2);
       C_obs_NE(i,j) = K0_NE./(1 + (d/a_NE).^2);
   end
   d = distance(lat0, long0, lat(i), long(i)) * 111 ; % km
   C_new_NN(i,1) = K0_NN./(1 + (d/a_NN).^2);
   C_new_EE(i,1) = K0_EE./(1 + (d/a_EE).^2);
   C_new_NE(i,1) = K0_NE./(1 + (d/a_NE).^2);
end

%% no correlation
if ismember('no corr', flags)
    C_obs_NE = zeros(p,p);
    C_new_NE = zeros(p,1);
end

if ismember('corr only', flags)
    C_obs_NN = zeros(p,p);
    C_obs_EE = zeros(p,p);
    C_new_NN = zeros(p,1);
    C_new_EE = zeros(p,1);
end

%% solve LSC
% Covariances
C_obs = [C_obs_NN, C_obs_NE; ...
         C_obs_NE, C_obs_EE];

C_new = [C_new_NN, C_new_NE; ...
         C_new_NE, C_new_EE];
     
% Observations  
V_obs = [Vn; Ve];

% Solve LSC
V_pred = C_new' * C_obs^-1 * V_obs;

end