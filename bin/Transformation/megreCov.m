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

for iiSite = 1:length(SplittedSites_list)
    iSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iiSite}))
        i_add = find(strcmp(names_SNX,SplittedSites.(SplittedSites_list{iiSite}){i}));
        if ~isempty(i_add)
            iSet = [iSet; i_add];
        end
    end
    if ~isempty(iSet) 
        %%
        for j=1:length(iSet)
            CovAveStack = zeros(3,3,j);
%             names_unified{iiSite} = SplittedSites_list{iiSite};
            for iiSite = [1:(iSet(1)-1),(iSet(end)+1):p]
                % avarage covariances b/w sonstrained sites w.r.t others, by columns   
                Cov2(iiSite, iSet(j)*3+0) = mean( Cov1(iiSite, iSet(j)*3+0) );
                Cov2(iiSite, iSet(j)*3+1) = mean( Cov1(iiSite, iSet(j)*3+1) );
                Cov2(iiSite, iSet(j)*3+2) = mean( Cov1(iiSite, iSet(j)*3+2) );
                Cov2(iSet(j)*3+0, iiSite) = mean( Cov1(iSet(j)*3+0, iiSite) );
                Cov2(iSet(j)*3+1, iiSite) = mean( Cov1(iSet(j)*3+1, iiSite) );
                Cov2(iSet(j)*3+2, iiSite) = mean( Cov1(iSet(j)*3+2, iiSite) );
                
                % average by block
                CovAveStack(1:3,1:3,j) = Cov1(iSet(j)*3-1+(1:3), iSet(j)*3-1+(1:3));
                CovMean = mean(CovAveStack,3);
                Cov2(iSet(j)*3-1+(1:3),iSet(j)*3-1+(1:3)) = CovMean;
            end
        end
        %%
        
    end
end

%% Simple take 1st only 
names_unified = names_SNX;


for iiSite = 1:length(SplittedSites_list)
    iSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iiSite}))
        i_add = find(strcmp(names_SNX,SplittedSites.(SplittedSites_list{iiSite}){i}));
        if ~isempty(i_add)
            iSet = [iSet; i_add];
        end
    end
%     iSet = nanclip(iSet)
    if ~isempty(iSet) 
        for j=1:length(iSet)
            iiSite = iSet(j);
            CRD_unified(iiSite,:) = mean(CRD(iSet,:),1);
            VEL_unified(iiSite,:) = mean(VEL(iSet,:),1);
            names_unified{iiSite} = SplittedSites_list{iiSite};
        end
    end
end

%%  megre stations with same 4-char name
uniq_names = unique(names_unified);
CRD_ave = zeros(length(uniq_names),3);
VEL_ave = zeros(length(uniq_names),3);

for i = 1:length(uniq_names)
    iSet = find(strcmp(names_unified,uniq_names(i)));
    CRD_ave(i,:) = mean(CRD_unified(iSet,:),1);
    VEL_ave(i,:) = mean(VEL_unified(iSet,:),1);
end



end