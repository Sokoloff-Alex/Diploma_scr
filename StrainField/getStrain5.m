function [Strain] = getStrain5(PointA, PointB, PointC, VelA, VelB, VelC)
% compute strain field
% Normal strain and Shear Strain b/w 2 points
%
% Alexandr Sokolov, KEG
% 28.11.2016

% block :   Bx 
%           AC    


vxAB = VelB(1) - VelA(1);
vyAB = VelB(2) - VelA(2);
vxAC = VelC(1) - VelA(1);
vyAC = VelC(2) - VelA(2);

xAB = deg2km( PointB(1) - PointA(1) )*1000; % [[m]
yAB = deg2km( PointB(2) - PointA(2) )*1000; %  [m]
xAC = deg2km( PointC(1) - PointA(1), 6378*cosd(PointA(2)) ) * 1000; % [m]
yAC = deg2km( PointC(2) - PointA(2), 6378*cosd(PointA(2)) ) * 1000; % [m]


%% Velocity gradient

% velocity
V = [vxAB, vyAB, vxAC, vyAC]';

% Desing matrix
A = [xAB yAB  0   0  
      0   0  xAB yAB
     xAC yAC  0   0
      0   0  xAC yAC];

% Velocity gradient Estimation
F = (A'*A)^-1*A'*V;

dvxdx = F(1);
dvxdy = F(2); 
dvydx = F(3);
dvydy = F(4);

dV = [     dvxdx,         1/2*(dvxdy+dvydx);
      1/2*(dvxdy+dvydx)        dvydy      ];
  
E = 1/2 * (dV + dV');
W = 1/2 * (dV - dV');

exx = E(1,1);
eyy = E(2,2);
exy = E(1,2);

w = 1/2 * (dvxdy - dvydx)

%% Alternative, Shorter way :

A2 = [xAB yAB  0   yAB;
       0  xAB yAB -xAB;
      xAC yAC  0   yAC;
       0  yAC yAC -xAC];
   
S = (A2'*A2)^-1*A2'*V;
S = [exx, exy, eyy, w];
    
%% Principal strain:

n1 = (exx+eyy)/2 + sqrt( ( (exx-eyy)/2 )^2 + exy^2 );
n2 = (exx+eyy)/2 - sqrt( ( (exx-eyy)/2 )^2 + exy^2 );

Theta_n1 = 1/2 * atan2d(2*exy , (exx-eyy));

NormalStrain = [n1, n2, Theta_n1];

%% Shear Strain

e12_max =  sqrt( ( (exx-eyy)/2 )^2 + exy^2 );

Theta_s = 1/2 * atan2d((eyy-exx) , (2*exy));

ShearStrain = [e12_max, Theta_s];

%% output result

x = (PointA(1) + PointC(1)) / 2;
y = (PointA(2) + PointB(2)) / 2;

Strain = [x, y, NormalStrain, ShearStrain, w];



end