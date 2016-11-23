function [StrainGrid] = getStrainMap2(Deformation)
% compute strain field
% result is shecker-board grid
%
% Alexandr Sokolov, KEG
% 14.11.2016

%%
LongStack = Deformation(:,1);
LatStack  = Deformation(:,2);
Vel       = Deformation(:,[3,4]);

% wrap stack into 2D grid
[VelGrid, LongGrid, LatGrid] =  vector2grid(Vel, LongStack, LatStack)

StepLong = LongGrid(1,2) - LongGrid(1,1)
StepLat  = LatGrid(2,1)  - LatGrid(1,1) 

sizeGrid = size(LongGrid);
StrainGrid = NaN([sizeGrid-1, 8]);

% points order  :  Bx
%                  AC


for iLong = 1:size(VelGrid,1)-1
   for iLat = 1:size(VelGrid,2)-1
       PointA = [LongGrid(iLat, iLong),     LatGrid(iLat,     iLong)];
       PointB = [LongGrid(iLat, iLong),     LatGrid(iLat + 1, iLong)];
       PointC = [LongGrid(iLat, iLong + 1), LatGrid(iLat,     iLong)];
       VelA = squeeze(VelGrid(iLat,    iLong,:)); 
       VelB = squeeze(VelGrid(iLat + 1,iLong,:));
       VelC = squeeze(VelGrid(iLat,    iLong + 1,:));
       
%        % check points order
%        [PointB(1), 0; ...
%         PointA(1), PointC(1) ]
%     
%     
%        [PointB(2), 0; ...
%         PointA(2), PointC(2) ]
    
       % get Strain b/w grid ponts
       StrainBlock = getStrain3(PointA, PointB, PointC, VelA, VelB, VelC);
       
       % put into grid
       StrainGrid(iLong, iLat, 1:8) = StrainBlock;
%        StrainGrid(iLong, iLat, 7:8) = [LongGrid(iLong, iLat) + StepLong/2, LatGrid(iLong, iLat) + StepLat/2];
   end
end   

end