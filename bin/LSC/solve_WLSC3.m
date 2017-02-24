function [V_pred, rmsFitting, V_noise_pred, Csig0] = solve_WLSC3(lat0, long0, lat, long, Venu, CovVel, Cov_scale, varargin)
% solve Least Square Collocation in 3D
%
% input   :     lat0, long0  - coordinates of grid point , [deg]
%               lat, long    - coordinates of observation points
%               Venu         - velocity components
%               CovVel       - covariance matrix of velosities from SNX COVA
%                              and transfrmed by covXYZ2ENU()
%               flags        - approx. function type ('exp1', 'Hirvonen' or 'exp2')
%               varargin     - 'plot' to plot covariance functions
%                            - '-v' verbose mode           
% output  :     V_pred       - interpolated / predicted velocity vector
%               V_noise_pred - predicted noise, propagated from measurement
%                              noise
%               rmsFitting   - rms of fitting Data points to to empirical Covariance curve
%
% Example :
%   
% V_pred = solve_LSC(iLat, iLong, lat, long, Ve, Vn, Vu,CovVel, flags, 'exp1','-v')
%
% Alexandr Sokolov, KEG
% 07.12.2016

flags = varargin{1};

% parse fType
fType = flags(ismember(flags, {'exp1', 'Hirvonen','normal','exp2','Gaussian','gaussian'}));

%
p = length(lat);

Ve = Venu(:,1);
Vn = Venu(:,2); 
Vu = Venu(:,3);

%% COMPUTE CORRELATIONS
[Corr_EE,Corr_NN,Corr_UU,Corr_EN,Corr_EU, Corr_NU, nClasses, scale] = getEmpiricalCorrelations3(lat, long, Ve, Vn, Vu);

%% get mean and trim columns with all NaN
[x_EE, y_EE] = mean_No_NaN(Corr_EE);
[x_NN, y_NN] = mean_No_NaN(Corr_NN);
[x_UU, y_UU] = mean_No_NaN(Corr_UU);
[x_EN, y_EN] = mean_No_NaN(Corr_EN);
[x_EU, y_EU] = mean_No_NaN(Corr_EU);
[x_NU, y_NU] = mean_No_NaN(Corr_NU);

%% Filter Noise 
% Cs(0) = Cl(0) - Cr(0)
if max(ismember(flags, {'filter', '-f'}))
    CrEE = trace(CovVel(0*p+1:1*p, 0*p+1:1*p)) * 1000^2 * 1 / p;
    CrNN = trace(CovVel(1*p+1:2*p, 1*p+1:2*p)) * 1000^2 * 1 / p;
    CrUU = trace(CovVel(2*p+1:3*p, 2*p+1:3*p)) * 1000^2 * 1 / p;
    
    Cr_TH = 0.05;
    if CrEE < Cr_TH && CrEE > 0.05*y_EE(1) ; CrEE = Cr_TH; end
    if CrNN < Cr_TH && CrNN > 0.05*y_NN(1) ; CrNN = Cr_TH; end
    if CrUU < Cr_TH && CrUU > 0.05*y_UU(1) ; CrUU = Cr_TH; end
    
    if CrEE < Cr_TH && CrEE < 0.05*y_EE(1) ; CrEE = 0.05*y_EE(1); end
    if CrNN < Cr_TH && CrNN < 0.05*y_NN(1) ; CrNN = 0.05*y_NN(1); end
    if CrUU < Cr_TH && CrUU < 0.05*y_UU(1) ; CrUU = 0.05*y_UU(1); end
    
    y_EE(1) = y_EE(1) - CrEE;
    y_NN(1) = y_NN(1) - CrNN;
    y_UU(1) = y_UU(1) - CrUU;  
end

%%
x_EE(2:end) = x_EE(2:end)*scale - scale/2;
x_NN(2:end) = x_NN(2:end)*scale - scale/2;
x_UU(2:end) = x_UU(2:end)*scale - scale/2;
x_EN(2:end) = x_EN(2:end)*scale - scale/2;
x_EU(2:end) = x_EU(2:end)*scale - scale/2;
x_NU(2:end) = x_NU(2:end)*scale - scale/2;


if y_EN(1) < 0
    b_EN =  -abs(max(y_EN));
else
    b_EN =   abs(min(y_EN));
end

if y_EU(1) < 0
    b_EU =  -abs(max(y_EU));
else
    b_EU =   abs(min(y_EU));
end

if y_NU(1) < 0
    b_NU =  -abs(max(y_NU));
else
    b_NU =   abs(min(y_NU));
end

%% add values to tail for better fitting
if ismember('tail', flags)
    n = 5;
    x_EE = [x_EE, x_EE(end)+scale:scale:x_EE(end)+scale*n];
    x_NN = [x_NN, x_NN(end)+scale:scale:x_NN(end)+scale*n];
    x_UU = [x_UU, x_UU(end)+scale:scale:x_UU(end)+scale*n];
    x_EN = [x_EN, x_EN(end)+scale:scale:x_EN(end)+scale*n];
    x_EU = [x_EU, x_EU(end)+scale:scale:x_EU(end)+scale*n];
    x_NU = [x_NU, x_NU(end)+scale:scale:x_NU(end)+scale*n];
    
    y_EE = [y_EE, min(y_EE)*ones(1,n)];
    y_NN = [y_NN, min(y_NN)*ones(1,n)];
    y_UU = [y_UU, min(y_UU)*ones(1,n)];
    y_EN = [y_EN,     -b_EN*ones(1,n)];
    y_EU = [y_EU,     -b_EU*ones(1,n)];
    y_NU = [y_NU,     -b_NU*ones(1,n)];
end

if ismember('tail 0', flags)
    n = 5;
    x_EE = [x_EE, x_EE(end)+scale:scale:x_EE(end)+scale*n];
    x_NN = [x_NN, x_NN(end)+scale:scale:x_NN(end)+scale*n];
    x_UU = [x_UU, x_UU(end)+scale:scale:x_UU(end)+scale*n];
    x_EN = [x_EN, x_EN(end)+scale:scale:x_EN(end)+scale*n];
    x_EU = [x_EU, x_EU(end)+scale:scale:x_EU(end)+scale*n];
    x_NU = [x_NU, x_NU(end)+scale:scale:x_NU(end)+scale*n];
    
    y_EE = [y_EE, zeros(1,n)];
    y_NN = [y_NN, zeros(1,n)];
    y_UU = [y_UU, zeros(1,n)];
    y_EN = [y_EN, zeros(1,n)];
    y_EU = [y_EU, zeros(1,n)];
    y_NU = [y_NU, zeros(1,n)];
end

%% fit curves
if ismember('bias', flags)
    if min(y_EE) < 0
        [coeff_EE, rmsEE] = fitCovar(fType, x_EE, y_EE + abs(min(y_EE)), [y_EE(1),-0.001]);
    else
        [coeff_EE, rmsEE] = fitCovar(fType, x_EE, y_EE, [y_EE(1), -0.001]);
    end
    
    if min(y_NN) < 0
        [coeff_NN, rmsNN] = fitCovar(fType, x_NN, y_NN + abs(min(y_NN)), [y_NN(1),-0.001]);
    else
        [coeff_NN, rmsNN] = fitCovar(fType, x_NN, y_NN, [y_NN(1), -0.001]);
    end
    
    if min(y_UU) < 0
        [coeff_UU, rmsUU] = fitCovar(fType, x_UU, y_UU + abs(min(y_UU)), [y_UU(1),-0.001]);
    else
        [coeff_UU, rmsUU] = fitCovar(fType, x_UU, y_UU, [y_UU(1), -0.001]);
    end
    
    
    [coeff_EN, rmsEN] = fitCovar(fType, x_EN, y_EN + b_EN,           [y_EN(1),-0.001]);
    [coeff_EU, rmsEU] = fitCovar(fType, x_EU, y_EU + b_EU,           [y_EU(1),-0.001]);
    [coeff_NU, rmsNU] = fitCovar(fType, x_NU, y_NU + b_NU,           [y_NU(1),-0.001]);
    
else
    [coeff_EE, rmsEE] = fitCovar(fType, x_EE, y_EE, [y_EE(1), -0.001]);
    [coeff_NN, rmsNN] = fitCovar(fType, x_NN, y_NN, [y_NN(1), -0.001]);
    [coeff_UU, rmsUU] = fitCovar(fType, x_UU, y_UU, [y_UU(1), -0.001]);
    [coeff_EN, rmsEN] = fitCovar(fType, x_EN, y_EN, [y_EN(1), -0.001]);
    [coeff_EU, rmsEU] = fitCovar(fType, x_EU, y_EU, [y_EU(1), -0.001]);
    [coeff_NU, rmsNU] = fitCovar(fType, x_NU, y_NU, [y_NU(1), -0.001]);
end

coeff_EE(1) = y_EE(1);
coeff_NN(1) = y_NN(1);
coeff_UU(1) = y_UU(1);
coeff_EN(1) = y_EN(1);
coeff_EU(1) = y_EU(1);
coeff_NU(1) = y_NU(1);

rmsFitting = [rmsEE, rmsNN, rmsUU, rmsEN, rmsEU, rmsNU];

%% optimize coeff b
if ismember('opt b', flags)
    b_min = 0.005;
    if (abs(coeff_EE(2)) > b_min)
        coeff_EE(2) = -b_min;
    end
    if (abs(coeff_NN(2)) > b_min)
        coeff_NN(2) = -b_min;
    end
    if (abs(coeff_UU(2)) > b_min)
        coeff_UU(2) = -b_min;
    end
    
    
end





%% plot Covariance functions, y=a*exp(b*d)
d = (0:(max(x_NN)/100):max(x_NN));
clr = lines(8);
if ismember('plot', flags)
    figure;
    hold on
    grid on
    title( ['point :: lat ', num2str(lat0), ...
        ' ; long  ',   num2str(long0),      ...
        ' ; # obs.: ', num2str(length(lat)),...
        ' ; step = ',  num2str(scale),      ...
        ' km', ';  # Classes ', num2str(nClasses)] );
    pl1 = plot(x_EE, y_EE, 'o--', 'Color', clr(1,:));
    pl2 = plot(x_NN, y_NN, 'o--', 'Color', clr(2,:));
    pl3 = plot(x_UU, y_UU, 'o--', 'Color', clr(3,:));
    pl4 = plot(x_EN, y_EN, 'o--', 'Color', clr(4,:)); 
    pl5 = plot(x_EU, y_EU, 'o--', 'Color', clr(5,:)); 
    pl6 = plot(x_NU, y_NU, 'o--', 'Color', clr(6,:)); 

    pl7  = plot(d, empiricalCovariance(fType, coeff_EE, d), 'Color', clr(1,:), 'LineWidth',2);
    pl8  = plot(d, empiricalCovariance(fType, coeff_NN, d), 'Color', clr(2,:), 'LineWidth',2);
    pl9  = plot(d, empiricalCovariance(fType, coeff_UU, d), 'Color', clr(3,:), 'LineWidth',2);
    pl10 = plot(d, empiricalCovariance(fType, coeff_EN, d), 'Color', clr(4,:), 'LineWidth',2);
    pl11 = plot(d, empiricalCovariance(fType, coeff_EU, d), 'Color', clr(5,:), 'LineWidth',2);
    pl12 = plot(d, empiricalCovariance(fType, coeff_NU, d), 'Color', clr(6,:), 'LineWidth',2);
    
    pl13 = plot(d, coeff_NN(1).*exp(-0.005*d ), '.-k');
    pl14 = plot(d, coeff_NN(1).*exp(-0.05*d),   '--k');

    pl15 = plot(0,y_EE(1) + CrEE, '*',  'Color', clr(1,:));
    pl16 = plot(0,y_NN(1) + CrNN, '*',  'Color', clr(2,:));
    pl17 = plot(0,y_UU(1) + CrUU, '*',  'Color', clr(3,:));

    
    legend([pl7 pl8 pl9 pl10 pl11 pl12 pl13 pl14], ...
        ['C_E_E C0 = ',  num2str(coeff_EE(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_EE(2),'%.3g'), ...
        ' ; rms_E_E = ', num2str(rmsEE)],...
        ['C_N_N C0 = ',  num2str(coeff_NN(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_NN(2),'%.3g'), ...
        ' ; rms_N_N = ', num2str(rmsNN)],...
        ['C_U_U C0 = ',  num2str(coeff_UU(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_UU(2),'%.3g'), ...
        ' ; rms_U_U = ', num2str(rmsUU)],...
        ['C_E_N  a = ',  num2str(coeff_EN(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_EN(2),'%.3g'), ...
        ' ; rms_N_N = ', num2str(rmsEN)],...
        ['C_E_U  a = ',  num2str(coeff_EU(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_EU(2),'%.3g'), ...
        ' ; rms_E_U = ', num2str(rmsEU)],...
        ['C_N_U  a = ',  num2str(coeff_NU(1),'%.2g'),  ...
               ' b = ',  num2str(coeff_NU(2),'%.3g'), ...
        ' ; rms_N_U = ', num2str(rmsNU)],...
        ['stable: C(d) = ' num2str(coeff_NN(1),'%.2g'),  ...
               '*exp(',num2str(-0.005,'%.2g'),'*d)'],  ...
        ['deform: C(d) = ', num2str(coeff_NN(1),'%.2g'),  ...
               '*exp(', num2str(-0.05,'%.2g'),'*d)']);
     xlabel('Classes, (distance [km])')
     ylabel('Covariance')
end

%% histogram of distances
if ismember('hist', flags)
    arc = zeros(p);
    dist = [];
    for i = 1:p
        arc(:,i) = distance(lat(i), long(i), lat, long) * 111 ; % km
        dist = [dist; arc(i:p,i)];
    end
    nelements = hist(ceil(dist/scale)*scale,length(x_NN));
    b2 = bar(x_NN, nelements/max(nelements)*max([y_NN(1), y_EE(1), y_NE(1)]),0.4);
    set(get(b2,'Children'),'FaceAlpha',0.2)
    text(x_NN, nelements/max(nelements)*max([y_NN(1), y_EE(1), y_NE(1)]),num2str(nelements'))
    xlabel('Classes , [km]')
    hold off
end

%% Build Covariance matrices
C_obs_EE = zeros(p,p);
C_obs_NN = zeros(p,p);
C_obs_UU = zeros(p,p);
C_obs_EN = zeros(p,p);
C_obs_EU = zeros(p,p);
C_obs_NU = zeros(p,p);

C_new_EE = zeros(p,1);
C_new_NN = zeros(p,1);
C_new_UU = zeros(p,1);
C_new_EN = zeros(p,1);
C_new_EU = zeros(p,1);
C_new_NU = zeros(p,1);


for i = 1:p
   for j = 1:p
       d = greatcircleArc(lat(i), long(i), lat(j), long(j)) * 111 ; % km
       C_obs_EE(i,j) = empiricalCovariance(fType, coeff_EE, d);
       C_obs_NN(i,j) = empiricalCovariance(fType, coeff_NN, d);
       C_obs_UU(i,j) = empiricalCovariance(fType, coeff_UU, d);
       C_obs_EN(i,j) = empiricalCovariance(fType, coeff_EN, d);
       C_obs_EU(i,j) = empiricalCovariance(fType, coeff_EU, d);
       C_obs_NU(i,j) = empiricalCovariance(fType, coeff_NU, d);
   end
   d = greatcircleArc(lat0, long0, lat(i), long(i)) * 111 ; % km  
   C_new_EE(i,1) = empiricalCovariance(fType, coeff_EE, d);
   C_new_NN(i,1) = empiricalCovariance(fType, coeff_NN, d);
   C_new_UU(i,1) = empiricalCovariance(fType, coeff_UU, d);
   C_new_EN(i,1) = empiricalCovariance(fType, coeff_EN, d);
   C_new_EU(i,1) = empiricalCovariance(fType, coeff_EU, d);
   C_new_NU(i,1) = empiricalCovariance(fType, coeff_NU, d);
end

%% no correlation
if ismember('no corr', flags)
    C_obs_EN = zeros(p,p);
    C_new_EN = zeros(p,1);
    C_obs_EU = zeros(p,p);
    C_new_EU = zeros(p,1);
    C_obs_NU = zeros(p,p);
    C_new_NU = zeros(p,1);
end

if ismember('corr only', flags)
    C_obs_EE = zeros(p,p);
    C_obs_NN = zeros(p,p);
    C_obs_UU = zeros(p,p);
    C_new_EE = zeros(p,1);
    C_new_NN = zeros(p,1);
    C_new_UU = zeros(p,1);
end

%% Covariances
C_obs = [C_obs_EE, C_obs_EN, C_obs_EU 
         C_obs_EN, C_obs_NN, C_obs_NU 
         C_obs_EU, C_obs_NU, C_obs_UU];

C_new = [C_new_EE, C_new_EN, C_new_EU
         C_new_EN, C_new_NN, C_new_NU
         C_new_EU, C_new_NU, C_new_UU];

%% Rescale Covariance Matrix accoring to SINEX
%  VARIANCE FACTOR                     1.800168208557666
%  [m/yr]^2  -> [mm/yr]^2              1000^2
%  Bernese scale                       10-100
Cnoise = CovVel;
% Cnoise([1:p],[1:p]) = Cnoise([1:p],[1:p])*1000^2;
Cnoise = Cnoise*1000^2 * Cov_scale;

%% Solve LSC
% Observations  
V_obs = [Ve - mean(Ve); Vn - mean(Vn); Vu - mean(Vu)];

V_pred = C_new' * (C_obs + Cnoise)^-1 * V_obs;

V_pred = ( V_pred + [mean(Ve); mean(Vn); mean(Vu)] )' ;

%% propagate noise
% SigmaSquare = diag(Cnoise);
% V_noise_pred = C_new' * (C_obs + Cnoise)^-1 * SigmaSquare;
% V_noise_pred = V_noise_pred';


%% Khale, error in dimentions
% V_noise_pred = C_new' * (C_obs + Cnoise)^-1 * C_new; %  * SigmaSquare; 
% V_noise_pred = diag(V_noise_pred)';

%% Mihkail E.
Cs0 = [coeff_EE(1),0,0; 
       0,coeff_NN(1),0;
       0,0,coeff_UU(1)];

V_noise_pred = Cs0 - C_new' * (C_obs + Cnoise)^-1 * C_new; 
V_noise_pred = sqrt(diag(V_noise_pred))';

Csig0 = sqrt(diag(Cs0));
% 
%% verbose mode
if ismember('-v', flags)
    disp(['point ::', ... 
          '; lat = ',  num2str(lat0, '%5.3f'), ...
          '; long = ', num2str(long0,'%5.3f'), ...
          ' # obs.: ', num2str(p), ...
          '; Ve = ', num2str(V_pred(1), '%#+7.4f'), ...
          '; Vn = ', num2str(V_pred(2), '%#+7.4f'), ...
          '; Vu = ', num2str(V_pred(3), '%#+7.4f'), ...
          '; SigmaVe^2 = ', num2str(V_noise_pred(1), '%#+7.4f'), ...
          '; SigmaVn^2 = ', num2str(V_noise_pred(2), '%#+7.4f'), ...
          '; SigmaVu^2 = ', num2str(V_noise_pred(3), '%#+7.4f'), ...
          ' ; rmsEE = ', num2str(rmsEE), ...
          ' ; rmsNN = ', num2str(rmsNN), ...
          ' ; rmsUU = ', num2str(rmsUU), ...
          ' ; rmsEN = ', num2str(rmsEN), ...
          ' ; rmsEU = ', num2str(rmsEU), ...
          ' ; rmsNU = ', num2str(rmsNU), ...
          ' :: Cee: b = ', num2str(coeff_EE(2), '%#7.4f'), ...
          ' :: Cnn: b = ', num2str(coeff_NN(2), '%#7.4f'), ...
          ' :: Cuu: b = ', num2str(coeff_UU(2), '%#7.4f'), ...
          ' :: Cen: b = ', num2str(coeff_EN(2), '%#7.4f'), ...
          ' :: Ceu: b = ', num2str(coeff_EU(2), '%#7.4f'), ...
          ' :: Cnu: b = ', num2str(coeff_NU(2), '%#7.4f') ]);
end

end