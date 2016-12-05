function Cov2 = extractCovariance(Cov1, iiSetOfBlocks, varargin )
% function to extract only part of covariance matrix Cov1
% by selected row and colunms of station defined by  iiSetOfBlocks
% and output into Cov2
% dim - dimenttions defined as array of dimentions to be used,
% like % [1 2], [1 2 3]
% and also separate dimentions: [[NN][NE]
%                                [NE][EE]]; by usign flag 'split'

dim  = varargin{1};
flag = varargin{2};

if isempty(varargin)
   dim = [1 2 3]; % use all 3 dimentions
end

iiElem = [];
for i = 1:length(dim)
    iiElem = [iiElem, (iiSetOfBlocks-1)*3 + dim(i) ];
end

if ismember(flag, {'split'})
%     disp('split')
else
%     disp('do not split')
    iiElem = sort(iiElem);
end
Cov2 = Cov1(:,iiElem); 
Cov2 = Cov2(iiElem,:);

end