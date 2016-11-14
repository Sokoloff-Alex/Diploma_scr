function [MapGrid, LongGrid, LatGrid] =  vector2grid(MapStack, LongStack, LatStack)
% function to wrap vector stack(p,1) into matrix(n,m), n * m = p;
%
% Alexandr Sokolov, KEG
% 14.11.2016

LongLength = length(unique(sort(LongStack)));
LatLength  = length(unique(sort(LatStack)));
MapNCol = size(MapStack,2);

LongGrid = NaN(LongLength, LatLength);
LatGrid  = NaN(LongLength, LatLength);
MapGrid  = NaN(LongLength, LatLength,MapNCol);

p = 0;
for i = 1:LongLength
    for j= 1:LatLength
       p = p + 1;
       MapGrid(i,j,:)  = MapStack(p,:);
       LongGrid(i,j) = LongStack(p);
       LatGrid(i,j)  = LatStack(p);
    end
end
    
end