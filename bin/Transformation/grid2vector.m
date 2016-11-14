function [vector] = grid2vector(grid)
% function to stratch grid/matrix into vector/stack
% do colunm-wise
%
% Input : grid/matrix,  [n,m]
% Output: vector/stack, [p,1], p = n*m
%
% command:
% [vector] = grid2vector(grid)
%
% Alexandr Sokolov, KEG
% 14.11.2016

[n,m] = size(grid);

for i = 1:m
    vector(i*n-(n-1):i*n,1) = grid(:,i); 
end

end
