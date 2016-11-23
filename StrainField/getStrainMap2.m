function [Strain] = getStrainMap2(Deformation)
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

Strain = NaN(size(VelGrid,1)*2-1, size(VelGrid,2)*2-1, 8);
% Grid   = NaN(size(VelGrid,1)*2-1, size(VelGrid,2)*2-1, 2);

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
       Strain = getStrain3(PointA,VelB-VelA,VelC-VelA);
       
       % put into grid
       Strain(iLong*2-1, iLat*2, 1:6) = StrainAB;
       Strain(iLong*2-1, iLat*2, 7:8) = [LongGrid(iLat, iLong), LatGrid( iLat, iLong) + StepLat/2];
       
       Strain(iLong*2, iLat*2-1, 1:6) = StrainAC;
       Strain(iLong*2, iLat*2-1, 7:8) = [LongGrid(iLat, iLong) + StepLong/2, LatGrid( iLat, iLong)];
   end
end   

% %% todo:: add rigth and bottom egdes !
% %% add rigth egde
% iLong = size(VelGrid,1);
% for iLat = 1:size(VelGrid,1)-1
%     PointA = [LongGrid(iLong),     LatGrid(iLat)];
%     PointB = [LongGrid(iLong),     LatGrid(iLat + 1)];
%     VelA = squeeze(VelGrid(iLong,   iLat,:)); 
%     VelB = squeeze(VelGrid(iLong,   iLat+1,:));
%     
%     % get Strain b/w grid ponts
%     StrainAB = getStrain3(PointA,VelA, PointB,VelB);
%     Strain(iLong*2-1, iLat*2,   1:6) = StrainAB;
%     Strain(iLong*2-1, iLat*2, 7:8) = [LongGrid(iLat, iLong), LatGrid( iLat, iLong) + StepLat/2];
% end
% 
% % add bottom edge
% iLat = size(VelGrid,2);
% for iLong = 1:size(VelGrid,2)-1;
%     PointA = [LongGrid(iLong),     LatGrid(iLat)];
%     PointC = [LongGrid(iLong + 1), LatGrid(iLat)];
%     VelA = squeeze(VelGrid(iLong,   iLat,:)); 
%     VelC = squeeze(VelGrid(iLong+1, iLat,:));
%     
%     % get Strain b/w grid ponts
%     StrainAC = getStrain3(PointA,VelA, PointC,VelC);
%     Strain(iLong*2,   iLat*2-1, 1:6) = StrainAC;
%     Strain(iLong*2, iLat*2-1, 7:8) = [LongGrid(iLat, iLong) + StepLong/2, LatGrid( iLat, iLong)];
% end

end