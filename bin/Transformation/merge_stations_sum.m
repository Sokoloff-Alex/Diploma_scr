function [Value_ave,uniq_names] = merge_stations_sum(Value,names)
% megre artificial stations
% by computing sum of Value, Tobs for example 
%

SplittedSites = load_splitted_sites();
SplittedSites_list = fieldnames(SplittedSites);

%%
names_unified = names;
Value_unified = Value;
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
            Value_unified(row,:) = mean(Value(iSet,:),1);
            names_unified{row} = SplittedSites_list{iSite};
        end
    end
end

%%  megre stations with same 4-char name
uniq_names = unique(names_unified);
Value_ave = zeros(length(uniq_names),size(Value,2));

for i = 1:length(uniq_names)
    iSet = find(strcmp(names_unified,uniq_names(i)));
    Value_ave(i,:) = sum(Value_unified(iSet,:),1);
end

end