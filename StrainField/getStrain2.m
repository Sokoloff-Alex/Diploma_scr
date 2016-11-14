function Strain = getStrain2(point1, vel1, point2, vel2)
% compute strain field
%
% Alexandr Sokolov, KEG
% 14.11.2016


x1 = point1(1);
y1 = point1(2);
x2 = point2(1);
y2 = point2(2);

%
dx = ged2km(x2-x1, (y1+y2)/2);
dy = deg2km(y2-y1);
dvx = vel2(1) - vel1(1);
dvy = vel2(2) - vel1(2);

%% get Velocity gradient
F = [dvx/dx , dvx/dy; 
     dvy/dx , dvy/dy];

fxx = F(1,1);
fyy = F(2,2);
fxy = F(1,2);
fyx = F(2,1);
    
    
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

%%
Strain = [NormalStrain, ShearStrain];
  

end