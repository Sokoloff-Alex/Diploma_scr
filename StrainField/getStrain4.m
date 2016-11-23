function [Strain] = getStrain4(PointA, PointB, PointC, VelA, VelB, VelC)
% compute strain field
% Normal strain and Shear Strain b/w 2 points
%
% Alexandr Sokolov, KEG
% 23.11.2016

% block :   Bx 
%           AC    

velAB = VelB - VelA;
velAC = VelC - VelA;

normAB1 = PointB(2) - PointA(2);                     % [m]
normAB2 = norm([0, PointB(2) - PointA(2)] + velAB ); % [m]

normAC1 = PointC(1) - PointA(1);                     % [m]
normAC2 = norm([PointC(1) - PointA(1), 0] + velAC);  % [m]


% Velocity gradient

fxx = ( normAC2 - normAC1 ) / normAC1;
fyy = ( normAB2 - normAB1 ) / normAB1;

fxy = velAB(1) / ( normAB1 + velAB(2));
fyx = velAC(2) / ( normAC1 + velAC(1));

F = [fxx fxy ; fyx, fyy];
    
%%  Rotation matrix R, Antisymmetric part
%   R = 1/2 * (F - F')
w =  (fyx - fxy)/2;

%%  Strain Tensor E, Symmetric part
%   E = 1/2 * (F + F')
E =  [fxx   , fxy+w; 
      fyx-w,  fyy ];
exx = E(1,1);
eyy = E(2,2);
exy = E(1,2); % = eyx

%% Principal Normal Strain
% solve for: det([exx-n, exy; eyx, eyy-n]) == 0! ;
n1 = (exx + eyy)/2 + sqrt( 1/4*(exx-eyy)^2 + exy^2 );
n2 = (exx + eyy)/2 - sqrt( 1/4*(exx-eyy)^2 + exy^2 );

alpha_n1 = 1/2 * atand(2*exy/(exx-eyy)); % alpha_n2 = alpha_n1 + 90;
alpha_n1 = alpha_n1 + 90;

NormalStrain = [n1, n2, alpha_n1];

%% Principal Shear Strain
s1 = + sqrt( 1/4*(exx-eyy)^2 + exy^2 ); % Only interesting in positive!
s2 = - sqrt( 1/4*(exx-eyy)^2 + exy^2 );

alpha_s1 = 1/2 * atand( -(exx-eyy)/(2*exy)); % alpha_s2 = alpha_s1 + 90;

ShearStrain = [s1, s2, alpha_s1];

%%
x = (PointA(1) + PointC(1)) / 2;
y = (PointA(2) + PointB(2)) / 2;
Strain = [x, y, NormalStrain, ShearStrain];
  

end