function [Cov2] = megreCov(Cov1, names_SNX)
% funiction to merge covariances for constrained stations
%
% Alexandr Sokolov, KEG
% 02.12.2016

Cov2 = Cov1;

p = size(Cov1,1)/3;

[SplittedSites] = load_splitted_sites();
SplittedSites_list = fieldnames(SplittedSites);

%% Average by columns constrained sites

setToBeMerged = [];

for iSite = 1:length(SplittedSites_list)
    iiSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iSite}))
        i_add = find(strcmp(names_SNX,SplittedSites.(SplittedSites_list{iSite}){i}));
        if ~isempty(i_add)
            iiSet = [iiSet; i_add];
        end
    end
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
        iAffected = setdiff(1:p,iiSet);
        for iRow = [iAffected]        
            for iSet = 1:length(iiSet)
                CovAveStackij(:,:,iSet) = Cov1((iRow-1)*3+(1:3), (iiSet(iSet)-1)*3+(1:3));
            end
            CovAveStackij_mean = mean(CovAveStackij,3);
            Cov2((iRow-1)*3+(1:3),(iiSet(1)-1)*3+(1:3)) = CovAveStackij_mean;
            Cov2((iiSet(1)-1)*3+(1:3),(iRow-1)*3+(1:3)) = CovAveStackij_mean; % Symmetric
        end
        setToBeMerged = [setToBeMerged, iiSet(1)];
    else
        setToBeMerged = [setToBeMerged, iSite];
    end
end

%% Simple take 1st only 
names_unified = names_SNX;


for iSite = 1:length(SplittedSites_list)
    iiSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iSite}))
        i_add = find(strcmp(names_SNX,SplittedSites.(SplittedSites_list{iSite}){i}));
        if ~isempty(i_add)
            iiSet = [iiSet; i_add];
        end
    end
%     iSet = nanclip(iSet)
    if ~isempty(iiSet) 
        for iSet=1:length(iiSet)
            iSite = iiSet(iSet);
            CRD_unified(iSite,:) = mean(CRD(iiSet,:),1);
            VEL_unified(iSite,:) = mean(VEL(iiSet,:),1);
            names_unified{iSite} = SplittedSites_list{iSite};
        end
    end
end

%%  megre stations with same 4-char name
uniq_names = unique(names_unified);
CRD_ave = zeros(length(uniq_names),3);
VEL_ave = zeros(length(uniq_names),3);

for i = 1:length(uniq_names)
    iiSet = find(strcmp(names_unified,uniq_names(i)));
    CRD_ave(i,:) = mean(CRD_unified(iiSet,:),1);
    VEL_ave(i,:) = mean(VEL_unified(iiSet,:),1);
end



end