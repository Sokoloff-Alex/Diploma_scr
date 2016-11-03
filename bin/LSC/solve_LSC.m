function [V_pred] = solve_LSC(lat0, long0, lat, long, Vn, Ve, fType, varargin)
% solve Least Square Collocation
%
% input   :     lat0, long0 - coordinates of grid point , [deg]
%               lat, long   - coordinates of observation points
%               Ve Vn       - velocity components
%               fType       - approx. function type ('exp1', 'Hirvonen' or 'exp2')
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
[Corr_NN, Corr_EE, Corr_NE, nClasses, scale] = getEmpiricalCorrelations2(lat, long, Vn, Ve);

%% get mean and trim columns with all NaN
[x_NN, y_NN] = mean_No_NaN(Corr_NN);
[x_EE, y_EE] = mean_No_NaN(Corr_EE);
[x_NE, y_NE] = mean_No_NaN(Corr_NE);

# x(1) = 0;
x_NN(2:end) = x_NN(e:end)*scale - scale/2;
x_EE(2:end) = x_EE(2:end)*scale - scale/2;
x_NE(2:end) = x_NE(2:end)*scale - scale/2;

if y_NE(1) < 0
    b_NE =  -abs(max(y_NE));
else
    b_NE =   abs(min(y_NE));
end

%% add values to tail for better fitting
if ismember('tail', flags)
    n = 5;
    x_NN = [x_NN, x_NN(end)+scale:scale:x_NN(end)+scale*n];
    x_EE = [x_EE, x_EE(end)+scale:scale:x_EE(end)+scale*n];
    x_NE = [x_NE, x_NE(end)+scale:scale:x_NE(end)+scale*n];
    
    y_NN = [y_NN, min(y_NN)*ones(1,n)];
    y_EE = [y_EE, min(y_EE)*ones(1,n)];
    y_NE = [y_NE,     -b_NE*ones(1,n)];
end

if ismember('tail 0', flags)
    n = 5;
    x_NN = [x_NN, x_NN(end)+scale:scale:x_NN(end)+scale*n];
    x_EE = [x_EE, x_EE(end)+scale:scale:x_EE(end)+scale*n];
    x_NE = [x_NE, x_NE(end)+scale:scale:x_NE(end)+scale*n];
    
    y_NN = [y_NN, zeros(1,n)];
    y_EE = [y_EE, zeros(1,n)];
    y_NE = [y_NE, zeros(1,n)];
end

%% fit curves
if ismember('bias', flags)
    coeff_NN = fitCovar(fType, x_NN, y_NN + abs(min(y_NN)), [y_NN(1),-0.001]);
    coeff_EE = fitCovar(fType, x_EE, y_EE + abs(min(y_EE)), [y_EE(1),-0.001]);
    coeff_NE = fitCovar(fType, x_NE, y_NE + b_NE,           [y_NE(1),-0.001]);   
else
    coeff_NN = fitCovar(fType, x_NN, y_NN, [y_NN(1), -0.001]);
    coeff_EE = fitCovar(fType, x_EE, y_EE, [y_EE(1), -0.001]);
    coeff_NE = fitCovar(fType, x_NE, y_NE, [y_NE(1), -0.001]);
end
coeff_NN(1) = y_NN(1);
coeff_EE(1) = y_EE(1);
coeff_NE(1) = y_NE(1);

%% verbose mode
if ismember('-v', flags)
    disp(['point ::', ... 
          ' lat = ',num2str(lat0, '%5.2f'), ...
          '; long = ',  num2str(long0,'%5.2f'), ...
          ' # obs.: ', num2str(p), ...
          ' :: Cnn: b = ', num2str(coeff_NN(2), '%#7.4f'), ...
          ' :: Cee: b = ', num2str(coeff_EE(2), '%#7.4f'), ...
          ' :: Cne: b = ', num2str(coeff_NE(2), '%#7.4f') ]);
end

%% plot Covariance functions, y=a*exp(b*d)
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
    pl4 = plot(d, empiricalCovariance(fType, coeff_NN, d), '-b', 'LineWidth',2);
    pl5 = plot(d, empiricalCovariance(fType, coeff_EE, d), '-r', 'LineWidth',2);
    pl6 = plot(d, empiricalCovariance(fType, coeff_NE, d), '-g', 'LineWidth',2);
    
    pl7 = plot(d, coeff_NN(1).*exp(-0.005*d ), '.-k');
    pl8 = plot(d, coeff_NN(1).*exp(-0.05*d), '--k');

    legend([pl4 pl5 pl6 pl7 pl8],['C_N_N C0 = ', num2str(coeff_NN(1),'%.2g'),  ...
                                        ' b = ', num2str(coeff_NN(2),'%.2g')], ...
                                 ['C_E_E C0 = ', num2str(coeff_EE(1),'%.2g'),  ...
                                        ' b = ', num2str(coeff_EE(2),'%.2g')], ...
                                 ['C_N_E  a = ', num2str(coeff_NE(1),'%.2g'),  ...
                                        ' b = ', num2str(coeff_NE(2),'%.2g')], ...
                                 ['stable: C(d) = ' num2str(coeff_NN(1),'%.2g'),  ...
                                        '*exp(',num2str(-0.005,'%.2g'),'*d)'],  ...
                                 ['deform: C(d) = ', num2str(coeff_NN(1),'%.2g'),  ...
                                        '*exp(', num2str(-0.05,'%.2g'),'*d)']);
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
       C_obs_NN(i,j) = empiricalCovariance(fType, coeff_NN, d);
       C_obs_EE(i,j) = empiricalCovariance(fType, coeff_EE, d);
       C_obs_NE(i,j) = empiricalCovariance(fType, coeff_NE, d);
   end
   d = distance(lat0, long0, lat(i), long(i)) * 111 ; % km  
   C_new_NN(i,1) = empiricalCovariance(fType, coeff_NN, d);
   C_new_EE(i,1) = empiricalCovariance(fType, coeff_EE, d);
   C_new_NE(i,1) = empiricalCovariance(fType, coeff_NE, d);
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
