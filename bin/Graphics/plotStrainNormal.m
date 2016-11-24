function plotStrainNormal(Strain, sc)
% funiction to plot Strain field on 2D plot
%
% input : Strain [Lambda1, Lambga2, Omega1] - data matrix for all points,
% Omega in [deg]
%                   
%         Long, Lat       - coord-s, of points
%         sc              - scale of arrows   
% 
% output:
%
% command:
% [fig] = plotStrain(Strain, Long, Lat)
%
% Alexandr Sokolov, DGFI
% 07.11.2016

[Stack_dilat, Stack_compr] = getStrain2quiver(Strain, sc);

%% Plot Strain map, faster

quiver(Stack_dilat(:,1),Stack_dilat(:,2),Stack_dilat(:,3),Stack_dilat(:,4), 0, 'b') % dilatiation
quiver(Stack_compr(:,1),Stack_compr(:,2),Stack_compr(:,3),Stack_compr(:,4), 0, 'r') % compression


%% plot strain map
% 
% for i = 1:size(Strain,1)
%     
%     % first axis
%     dx_L1 = L1(i)*cosd(Omega(i));
%     dy_L1 = L1(i)*sind(Omega(i));
% 
%     % second axis
%     dx_L2 = L2(i)*cosd(Omega(i)+90);
%     dy_L2 = L2(i)*sind(Omega(i)+90);
% 
%    if L1(i) > 0 
%       quiver(Long(i), Lat(i),  dx_L1*sc,  dy_L1*sc, 0, 'b', 'lineWidth',2) 
%       quiver(Long(i), Lat(i), -dx_L1*sc, -dy_L1*sc, 0, 'b', 'lineWidth',2) 
%    else 
%       quiver(Long(i) + dx_L1*sc, Lat(i) + dy_L1*sc, -dx_L1*sc, -dy_L1*sc, 0, 'r', 'lineWidth',2)   
%       quiver(Long(i) - dx_L1*sc, Lat(i) - dy_L1*sc, +dx_L1*sc, +dy_L1*sc, 0, 'r', 'lineWidth',2)   
%    end 
% 
%    if L2(i) > 0 
%       quiver(Long(i), Lat(i),  dx_L2*sc,  dy_L2*sc, 0, 'b', 'lineWidth',2) 
%       quiver(Long(i), Lat(i), -dx_L2*sc, -dy_L2*sc, 0, 'b', 'lineWidth',2) 
%    else 
%       quiver(Long(i) + dx_L2*sc, Lat(i) + dy_L2*sc, -dx_L2*sc, -dy_L2*sc, 0, 'r', 'lineWidth',2)   
%       quiver(Long(i) - dx_L2*sc, Lat(i) - dy_L2*sc, +dx_L2*sc, +dy_L2*sc, 0, 'r', 'lineWidth',2)   
%    end 
% end

end