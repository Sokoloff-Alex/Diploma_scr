function [Strain] = getStrainMap(Deformation)
% compute strain field
% result is shecker-board grid
%
% Alexandr Sokolov, KEG
% 14.11.2016

%%
LongStack = Deformation(:,1);
LatStack  = Deformation(:,2);
Vel = Deformation(:,[3,4]);

% wrap stack into 2D grid
[VelGrid, LongGrid, LatGrid] =  vector2grid(Vel, LongStack, LatStack);

Strain = NaN(size(VelGrid,1)*2-2, size(VelGrid,2)*2-2, 6);
Grid   = NaN(size(VelGrid,1)*2-2, size(VelGrid,2)*2-2, 2);

StepLong = LongGrid(1,2) - LongGrid(1,1)
StepLat  = LatGrid(2,1)  - LatGrid(1,1) 

for iLong = 1:size(VelGrid,1)-1
   for iLat = 1:size(VelGrid,2)-1
       PointA = [LongGrid(iLong),     LatGrid(iLat)];
       PointB = [LongGrid(iLong),     LatGrid(iLat + 1)];
       PointC = [LongGrid(iLong + 1), LatGrid(iLat)];
       VelA = squeeze(VelGrid(iLong,   iLat,:)); 
       VelB = squeeze(VelGrid(iLong,   iLat+1,:));
       VelC = squeeze(VelGrid(iLong+1, iLat,:));
       
       % get Strain b/w grid ponts
       StrainAB = getStrain2(PointA,VelA, PointB,VelB);
       StrainAC = getStrain2(PointA,VelA, PointC,VelC);
       
       % put into grid
       Strain(iLong*2,  iLat*2-1,1:6) = StrainAC;
       Strain(iLong*2-1,  iLat*2,7:8) = [LongGrid(iLat, iLong) + StepLong/2, LatGrid( iLat, iLong)];

       Strain(iLong*2-1,iLat*2,  1:6) = StrainAB;
       Strain(iLong*2,  iLat*2-1,7:8) = [LongGrid(iLat, iLong), LatGrid( iLat, iLong) + StepLat/2];
   end
end   

%% todo:: add left and upper egdes !


end