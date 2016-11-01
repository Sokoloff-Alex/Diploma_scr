function [iSet1, iSet0] = selectRange(array1,array2)
% provide range of indexes (iSet) of cells
% from array1 that appears in array2
% datatape in cell arays must be a String
%
% Command:
% [iSet1, iSet0] = selectRange(array1,array2)
% Alexandr Sokolov, KEG
% 13.10.2016

    iSet = [];
    for i = 1:length(array2)                          
        i_add = find(strcmp(array1, array2(i)))';
        iSet = [iSet, i_add];
    end
    iSet1  = sort(iSet);
    
%% same
    r = 1:length(array1);
    iSet2 = r(ismember(array1, array2) == 1);
    iSet0 = r(ismember(array2, array2) == 0);

end
