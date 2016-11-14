function [Strain, Grid] = getStrainMap(Deformation)
% compute strain field
%
% Alexandr Sokolov, KEG
% 14.11.2016

%%
Vel = Deformation(:,[3,4]);
LongStack = Deformation(:,1);
LatStack  = Deformation(:,2);


[VelGrid, LongGrid, LatGrid] =  vector2grid(Vel, LongStack, LatStack);

Strain = NaN(size(VelGrid)*2);
Grid   = NaN(size(VelGrid)*2);


for iLong = 1:(size(VelGrid,1)-1)
   for iLat = 1:(size(VelGrid,2)-1)
       pointA = [LongGrid(iLong),     LatGrid(iLat)];
       pointB = [LongGrid(iLong),     LatGrid(iLat + 1)];
       pointC = [LongGrid(iLong + 1), LatGrid(iLat)];
       velA = squeeze(VelGrid(iLong,     iLat)); 
       velB = squeeze(VelGrid(iLong,     iLat + 1)); 
       velC = squeeze(VelGrid(iLong + 1, iLat)); 
       StrainAB = getStrain2(pointA,velA, pointB,velB);
       StrainAC = getStrain2(pointA,velA, pointC,velC);
       Strain(iLong*2,  iLat*2-1)= StrainAC;
       Strain(iLong*2-1,iLat*2)  = StrainAB;
       
       Grid(iLong*2,ILat*2-1,1) = LatGrid(iLat,iLong);       
       Grid(iLong*2-1,ILat*2,2) = LongGrid(iLat,iLong);       
    
   end
end
       

end