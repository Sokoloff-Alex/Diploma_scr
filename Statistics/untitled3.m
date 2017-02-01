

SINEX = readSNX('../STA/FMC_IGB_W7.SNX');
covA = SINEX.SOLUTION.COVA_APRIORI;
covE = SINEX.SOLUTION.COVA_ESTIMATE;

%%
close all
figure(2)
hold on
imagesc(flipud(covE))
% imagesc(flipud(covE([1120:1140],[1120:1140])))

%%

l = size(covE)

l/6

SINEX.SITE.ID.CODE([10,187:190],:)

