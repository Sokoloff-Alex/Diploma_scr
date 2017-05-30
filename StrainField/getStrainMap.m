function [StrainStack] = getStrainMap(Deformation)
% compute strain field over the entire grid
%
% % the input data must be taken from regular grid
%
% input: coordinates(long, lat in degrees for 3 points) and values of 
% horizontal velocity components in [m/s]
%
% ouput: [x, y, n1, n2, Theta_n1, e12_max, Theta_s, w]
% result is Stack [p,7], where p is a number of grid points.
%
% with columns: x        - is a longitude [deg];
%               y        - is a latitude [deg];
%               n1       - is a principal strain component in first dimenton
%               n2       - is a principal strain component in second dimenton, n2 aothogonal to n1;
%               Theta_n1 - is an azimut of n1, [deg]
%               e12_max  - is a shear strain component;
%               Theta_s  - is an azimut of shera strain [deg];
%               w        - is a wedge;
%
% Alexandr Sokolov, KEG
% 23.11.2016

%% wrap stack into 2D grid
Grid = stack2grid(Deformation);

LongGrid = squeeze(Grid(:,:,1));
LatGrid  = squeeze(Grid(:,:,2));
VelGrid =  squeeze(Grid(:,:,3:end));

nLong = size(Grid,2);
nLat  = size(Grid,1); 

StrainStack = NaN((nLat-1)*(nLong-1),8);

% points order  :  Bx
%                  AC

%% compute Strain

row = 0;
for iLong = 1:nLong-1
   for iLat = 1:nLat-1
       PointA = [LongGrid(iLat, iLong),     LatGrid(iLat,     iLong)];
       PointB = [LongGrid(iLat, iLong),     LatGrid(iLat + 1, iLong)];
       PointC = [LongGrid(iLat, iLong + 1), LatGrid(iLat,     iLong)];
       VelA = squeeze(VelGrid(iLat,    iLong,:)); 
       VelB = squeeze(VelGrid(iLat + 1,iLong,:));
       VelC = squeeze(VelGrid(iLat,    iLong + 1,:));
    
       % get Strain b/w grid ponts
       StrainBlock = getStrain(PointA, PointB, PointC, VelA, VelB, VelC);
       
       % Save result in Stack
       row = row + 1;
       StrainStack(row,:) = StrainBlock;
   end
end   


end