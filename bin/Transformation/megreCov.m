function [Cov2Red] = megreCov(Cov1, names_SNX)
% funiction to merge covariances for constrained stations
%
% Alexandr Sokolov, KEG
% 02.12.2016

Cov2 = Cov1;

p = size(Cov1,1)/3;

[SplittedSites] = load_splitted_sites();
SplittedSites_list = fieldnames(SplittedSites);

%% Average by columns constrained sites

setToBeSkipped = [];
names_unified = names_SNX;

for iSite = 1:length(SplittedSites_list)
    iiSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iSite}))
        i_add = find(strcmp(names_SNX,SplittedSites.(SplittedSites_list{iSite}){i}));
        if ~isempty(i_add)
            iiSet = [iiSet; i_add];
        end
    end
    
    % average blocks
    if ~isempty(iiSet) 
        % average blocks at diagonal
        CovAveStackDiag = zeros(3,3,length(iiSet));
        for iSet=1:length(iiSet)
            % names_unified{iiSite} = SplittedSites_list{iiSite};
            iBlock = (iiSet(iSet)-1)*3+(1:3);
            CovAveStackDiag(:,:,iSet) = Cov1(iBlock,iBlock);
        end
        CovMeanDiag = mean(CovAveStackDiag,3);
        Cov2((iiSet(1)-1)*3+(1:3),(iiSet(1)-1)*3+(1:3)) = CovMeanDiag;

        % avarage covariances b/w constrained sites w.r.t others, by columns   
        CovAveStackij = zeros(3,3,length(iiSet));
        iAffected = setdiff(1:p,iiSet); % affected other sites
        for iRow = [iAffected]        
            for iSet = 1:length(iiSet)
                CovAveStackij(:,:,iSet) = Cov1((iRow-1)*3+(1:3), (iiSet(iSet)-1)*3+(1:3));
                names_unified{iiSet(iSet)} = SplittedSites_list{iSite};
            end
            CovAveStackij_mean = mean(CovAveStackij,3);
            Cov2((iRow-1)*3+(1:3),(iiSet(1)-1)*3+(1:3)) = CovAveStackij_mean;
            Cov2((iiSet(1)-1)*3+(1:3),(iRow-1)*3+(1:3)) = CovAveStackij_mean; % Symmetric
        end
        setToBeSkipped = [setToBeSkipped; iiSet(2:end)];
    end
end
names_unified;
setToBeSkipped;


%%  megre stations with same 4-char name
uniq_names = unique(names_unified);
p = length(uniq_names);
Cov2_red = zeros(3*p,3*p);

% prepare indices to preserve

% todo : average comm 4 char
iiSetSave = [];
for i = 1:p
    iiSet = find(strcmp(names_unified,uniq_names(i)));
    iiSetSave = [iiSetSave; iiSet(1)];
end


iiElemSave = sort([(iiSetSave-1)*3+1,(iiSetSave-1)*3+2, (iiSetSave-1)*3+3]);

Cov2Red = Cov2(:,iiElemSave);
Cov2Red = Cov2Red(iiElemSave,:);



end