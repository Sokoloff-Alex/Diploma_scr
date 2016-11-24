function [GridMatrix] =  stack2grid(Stack)
% function to wrap vector stack(p,1) into matrix(n,m), n * m = p;
%
% Alexandr Sokolov, KEG
% 23.11.2016

LongLength = length(unique(sort( Stack(:,1) )));
LatLength  = length(unique(sort( Stack(:,2) )));
StackNCol  = size(Stack,2);

GridMatrix = NaN(LatLength, LongLength, StackNCol);

p = 0;
for iLong = 1:LongLength
    for iLat= 1:LatLength
       p = p + 1;
       GridMatrix(iLat, iLong, :) = Stack(p,:);
    end
end

end