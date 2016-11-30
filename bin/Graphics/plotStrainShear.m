function [pl1, pl2, pl3, pl4] = plotStrainShear(Strain, scale)
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
w =          Strain(:,8);

%% convert data to quiver plot
% main axis:
dx_S1 = ShearMax.*cosd(AlphaShear);
dy_S1 = ShearMax.*sind(AlphaShear);

% split: left-lateral (CCW), rigth-lateral(CW)

N = length(w);
n = 1:N;
setCCW = n(w >= 0);
setCW  = n(w <  0);





%% plot strain map
% left-lateral, (CCW)
pl1 = quiver(Long(setCCW), Lat(setCCW),  dx_S1(setCCW)*scale,  dy_S1(setCCW)*scale, 0, 'b', 'ShowArrowHead','off', 'Color', [0,100,0]/255, 'LineWidth',2);
pl2 = quiver(Long(setCCW), Lat(setCCW), -dx_S1(setCCW)*scale, -dy_S1(setCCW)*scale, 0, 'b', 'ShowArrowHead','off', 'Color', [0,100,0]/255, 'LineWidth',2);

% rigth-lateral, (CW)
pl3 = quiver(Long(setCW), Lat(setCW),  dx_S1(setCW)*scale,  dy_S1(setCW)*scale, 0, 'b', 'ShowArrowHead','off', 'Color', [153 50 255]/255, 'LineWidth',2);
pl4 = quiver(Long(setCW), Lat(setCW), -dx_S1(setCW)*scale, -dy_S1(setCW)*scale, 0, 'b', 'ShowArrowHead','off', 'Color', [153 50 255]/255, 'LineWidth',2);

end