function [CRD_ave,VEL_ave,uniq_names] = merge_stations_precision(CRD,VEL,names)
% megre artificial stations
% by computing average value CRD and VEL
% applicable, since velocities are tigthly constraind and CRD only for map.
% 
%

SplittedSites = load_splitted_sites();
SplittedSites_list = fieldnames(SplittedSites);

%%
names_unified = names;
CRD_unified = CRD;
VEL_unified = VEL;
clc
for iSite = 1:length(SplittedSites_list)
    iSet = [];
    for i = 1:length(SplittedSites.(SplittedSites_list{iSite}))
        i_add = find(strcmp(names,SplittedSites.(SplittedSites_list{iSite}){i}));
        if ~isempty(i_add)
            iSet = [iSet; i_add];
        end
    end
%     iSet = nanclip(iSet)
    if ~isempty(iSet) 
        for j=1:length(iSet)
            row = iSet(j);
            CRD_unified(row,:) = mean(CRD(iSet,:),1);
            sigmaVe_est = sqrt( 1/sum( VEL(iSet,1).^2).^-1 );
            sigmaVn_est = sqrt( 1/sum( VEL(iSet,2).^2).^-1 );
            sigmaVu_est = sqrt( 1/sum( VEL(iSet,3).^2).^-1 );
            VEL_unified(row,:) = [sigmaVe_est, sigmaVn_est, sigmaVu_est];
            names_unified{row} = SplittedSites_list{iSite};
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
    
    % Estimation theory (Bamler R.), Lec.: Parameter Estimation Basics, page 21
    sigmaVe_est = sqrt(1/sum( VEL_unified(iSet,1).^2).^-1 );
    sigmaVn_est = sqrt(1/sum( VEL_unified(iSet,2).^2).^-1 );
    sigmaVu_est = sqrt(1/sum( VEL_unified(iSet,3).^2).^-1 );
    VEL_ave(i,:) = [sigmaVe_est, sigmaVn_est, sigmaVu_est];
end

end