function [pl1, pl2] = plotStrainShear(Strain, scale)
% funiction to plot Strain Shear field on 2D plot
%
% input : Strain [Long, lat, Lambda1, Lambda1 Omega1, s1, Alpha_s1] - data matrix for all points,
%                 Long, Lat Alpha, Omega in [deg]
%
%         scale           - scale of arrows   
% 
% command:
%   plotStrain(Strain, Long, Lat)
%
% Alexandr Sokolov, DGFI
% 23.11.2016


Long       = Strain(:,1); 
Lat        = Strain(:,2);
ShearMax   = Strain(:,6);
AlphaShear = Strain(:,7);

%% convert data to quiver plot
% main axis:
dx_S1 = ShearMax.*cosd(AlphaShear);
dy_S1 = ShearMax.*sind(AlphaShear);


%% plot strain map
pl1 = quiver(Long, Lat,  dx_S1*scale,  dy_S1*scale, 0, 'b', 'ShowArrowHead','off', 'Color',[0,100/255,0]);
pl2 = quiver(Long, Lat, -dx_S1*scale, -dy_S1*scale, 0, 'b', 'ShowArrowHead','off', 'Color',[0,100/255,0]);

end