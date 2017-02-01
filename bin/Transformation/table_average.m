function [table_out] = table_average(table_in, fieldNbr)
% function to average values for same station

FieldNames = fieldnames(table_in);
NumberOfFields = length(FieldNames);

uniq_Entries_list = unique(sort(table_in.Site));

NumberOfUniqEnrties = length(uniq_Entries_list);

CRD = zeros(NumberOfUniqEnrties,3);
VEL = zeros(NumberOfUniqEnrties,3);


for i = 1:NumberOfUniqEnrties
   CRD(i,:) = mean(table_in.CRD(ismember(table_in.Site, uniq_Entries_list(i) ),: ), 1); 
   VEL(i,:) = mean(table_in.VEL(ismember(table_in.Site, uniq_Entries_list(i) ),: ), 1); 
end

Site = uniq_Entries_list;
table_out = table(Site, CRD, VEL,  'RowNames',Site);

end
