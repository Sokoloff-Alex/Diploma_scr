function [Strain] = getStrain3(PointA, PointB, PointC, VelA, VelB, VelC)
% compute strain field
% Normal strain and Shear Strain b/w 2 points
%
% Alexandr Sokolov, KEG
% 23.11.2016

% block :   Bx 
%           AC    


velAB = VelB - VelA;
velAC = VelC - VelA;

normAB1 =          deg2km( PointB(2) - PointA(2) ) * 1000;            % [m]
normAB2 = norm([0; deg2km( PointB(2) - PointA(2) ) * 1000] + velAB ); % [m]

normAC1 =       deg2km( PointC(1) - PointA(1), 6378*cosd(PointA(2)) ) * 1000;               % [m]
normAC2 = norm([deg2km( PointC(1) - PointA(1), 6378*cosd(PointA(2)) ) * 1000; 0] + velAC ); % [m]


%% Velocity gradient
fxx = ( normAC2 - normAC1 ) / normAC1;
fyy = ( normAB2 - normAB1 ) / normAB1;

fxy = velAB(1) / ( normAB1 + velAB(2));
fyx = velAC(2) / ( normAC1 + velAC(1));

% F = [fxx fxy ; fyx, fyy];
    
%%  Rotation matrix R, Antisymmetric part
%   R = 1/2 * (F - F')
w = (fyx - fxy)/2;

%%  Strain Tensor E, Symmetric part
%   E = 1/2 * (F + F')
E =  [fxx   , fxy+w; 
      fyx-w,  fyy ];
exx = E(1,1);
eyy = E(2,2);
exy = E(1,2); % = eyx
% eyx = E(2,1);
%% Principal Normal Strain
% solve for: det([exx-n, exy; eyx, eyy-n]) == 0! ;
n1 = (exx + eyy)/2 + sqrt( 1/4*(exx-eyy)^2 + exy^2 );
n2 = (exx + eyy)/2 - sqrt( 1/4*(exx-eyy)^2 + exy^2 );

% alpha_n1 = 1/2 * atan2d(2*exy, (exx-eyy)); % alpha_n2 = alpha_n1 + 90;

Omega1   = atan2d(-(exx - n1) , exy ); % angle for n1

% if alpha_n1 ~= Omega1
%     [alpha_n1, Omega1]
% end

NormalStrain = [n1, n2, Omega1];

%% Principal Shear Strain
s1 = + sqrt( 1/4*(exx-eyy)^2 + exy^2 );   % Only interesting in positive (Maximum)!
% s2 = - sqrt( 1/4*(exx-eyy)^2 + exy^2 ); % s1 = -s2

alpha_s1 = 1/2 * atand( -(exx-eyy) / (2*exy));

% alpha_s1 =  atan2d(  exy , (exx - s1));

ShearStrain = [s1, alpha_s1];

%% get into stack
x = (PointA(1) + PointC(1)) / 2;
y = (PointA(2) + PointB(2)) / 2;
Strain = [x, y, NormalStrain, ShearStrain, w];
  

end