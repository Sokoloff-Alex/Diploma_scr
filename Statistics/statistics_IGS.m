

Available  = zeros(16,366,25);

IGS    = {'AJAC','BRST','BZRG','FFMJ','GENO','GRAS','GRAZ','HFL2','HFLK','HUEG','IENG','LROC','MARS','MEDI','OBE2','OBE4','PADO','POTS','TLSE','VENE','WTZR','ZIM2','ZIMM'};

%%

for iSite = 9
    siteREF = IGS{iSite};


    for i = 1:length(datapoolIGStimesort)   
        site = datapoolIGStimesort{i,1};
        iDDD = datapoolIGStimesort{i,2};
        iYY  = datapoolIGStimesort{i,3};
        if strcmpi(site, siteREF)
            Available(iYY,iDDD,iSite) = 1;
        end 
    end


    for i = 1:length(processedlistIGS)   
        site = processedlistIGS{i,3};
        iDDD = processedlistIGS{i,2};
        iYY  = processedlistIGS{i,1};
        if strcmpi(site, siteREF)
            Available(iYY,iDDD,iSite) = 2;
        end 
    end

%%
    close all
    
    figure(1)
    hold on
    title([siteREF, ', green: available, Brown: processed'])
    imagesc(squeeze(Available(:,:,iSite)))
    xlabel('DOY')
    ylabel('Year')
    xlim([1 366])
    ylim([4 16])
    pause(1)
end


