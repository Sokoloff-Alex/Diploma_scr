function Strain = getStrain(Deformation)
% compute strain field
%
% Alexandr Sokolov, KEG
% 14.11.2016

%%
Vel = Deformation(:,[3,4]);
LongGrid = Deformation(:,1);
LatGrid   = Deformation(:,2);

%%
N = length(LatGrid);

x1 = LongGrid; % x
y1 = LatGrid;  % y
x2 = zeros(N,1);
y2 = zeros(N,1);

NormalStrain = NaN(N,3);
ShearStrain  = NaN(N,3);

for i = 1:N
    %% todo: to Correct
    y2 = y1(i) + km2deg(Vel(i,1)/1000);
    x2 = x1(i) + km2deg(Vel(i,2)/1000, 6371*cosd(y1(i)) );    

    dx = x2 - x1(i); % [deg]
    dy = y2 - y1(i); % [deg]

    % Deformation, F = [fxx, fxy ; fyx, fyy];
    fxx = dx / x1(i);
    fyy = dy / y1(i);
    fxy = dx / (y1(i) + dy);
    fyx = dy / (x1(i) + dx);
%% todo end
    
    %%  Rotation matrix R, Antisymmetric part
    % R = 1/2 * (F - F')
    w = (fyx - fxy)/2;

    %%  Strain Tensor E,  E = 1/2*(F + F')  = [exx exy; eyx, eyy]; Symmetric part
    E =  [fxx   , fxy+w; 
          fyx-w,  fyy ];
    exx = E(1,1);
    eyy = E(2,2);
    exy = E(1,2); % = eyx

    %% Principal Normal Strain
    e1 = (exx + eyy)/2 + sqrt( 1/4*(exx-eyy)^2 + exy^2 );
    e2 = (exx + eyy)/2 - sqrt( 1/4*(exx-eyy)^2 + exy^2 );
    
    alpha_e1 = 1/2 * atand(2*exy/(exx-eyy)); % alpha_e2 = alpha_e1 + 90;
    
    NormalStrain(i,:) = [e1, e2, alpha_e1];

    %% Principal Shear Strain
    s1 = + sqrt( 1/4*(exx-eyy)^2 + exy^2 );
    s2 = - sqrt( 1/4*(exx-eyy)^2 + exy^2 );
    
    alpha_s1 = 1/2 * atand( -(exx-eyy)/(2*exy)); % alpha_s2 = alpha_s1 + 90;
    
    ShearStrain(i,:) = [s1, s2, alpha_s1];
    
end

Grid = [LongGrid, LatGrid];
Strain = struct('NormalStrain', NormalStrain, 'ShearStrain',ShearStrain, 'Grid',Grid);

end