function [Stack] = grid2stack(Grid)
% function to stack grid/matrix into stack
% do colunm-wise
%
% Input : grid/matrix,  [n,m,x]
% Output: vector/stack, [p,x], p = n*m
%
% command:
% [stack] = grid2stack(grid)
%
% Alexandr Sokolov, KEG
% 23.11.2016

[n,m,x] = size(Grid);

for i = 1:m
    for k = 1:x
        Stack( i*n-(n-1):i*n,k) = squeeze(Grid(:,i,k));    
    end
end

end
