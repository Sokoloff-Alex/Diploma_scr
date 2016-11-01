function [CRD_ave,VEL_ave,uniq_names] = merge_stations(CRD,VEL,names)
% megre artificial stations
% by computing average value CRD and VEL
% applicable, since velocities are tigthly constraind and CRD only for map.
% 
%
%%

%% merge stations with names of type: AA1A, digit in 3rd charecter
% table is overloaded for completeness

AJAC = {{'AJ1C','AJ2C','AJ3C'}};
AUTN = {{'AU1N','AU2N','AU3N'}};
AXPV = {{'AX1V','AX2V','AX3V','AX4V'}};
BRST = {{'BR1T','BR2T','BR3T','BR4T','BR5T','BR6T'}};
BSCN = {{'BS1N','BS2N','BS3N'}};
BZRG = {{'BZ1G','BZ2G','BZ4G','BZ5G','BZ6G'}};
CHIZ = {{'CH1Z'}};
COMO = {{'CO1O','CO2O'}};
DRES = {{'DR1S','DR2S','DR3S','DR4S'}};
EGLT = {{'EG1T','EG2T','EG3T'}};
ENTZ = {{'EN1Z','EN2Z','EN3Z'}};
GARI = {{'GA1I','GA2I'}};
GENO = {{'GE1O'}};
GRAS = {{'GR1S','GR2S','GR3S'}};
GRAZ = {{'GR1Z','GR2Z','GR3Z','GR4Z','GR5Z'}};
GSR1 = {{'GS11','GS21','GS31','GS41'}};
IENG = {{'IE1G'}};
IGMI = {{'IG1I','IG2I','IG3I'}};
KARL = {{'KA1L','KA2L','KA3L','KA4L'}};
LINZ = {{'LI1Z','LI2Z','LI3Z'}};
LROC = {{'LR1C'}};
MARS = {{'MA1S','MA2S','MA3S','MA4S','MA5S'}};
MEDI = {{'ME1I','ME2I'}};
MLVL = {{'ML1L','ML2L','ML3L'}};
MOPS = {{'MO1S','MO2S','MO3S'}};
PFA2 = {{'PF12'}};
PORE = {{'PO1E'}};
POTS = {{'PO1S','PO2S','PO3S'}};
PRAT = {{'PR1T','PR2T'}};
PUYV = {{'PU1V','PU2V'}};
ROVE = {{'RO1E','RO2E','RO3E'}};
SBG2 = {{'SB12'}};
SJDV = {{'SJ1V','SJ2V'}};
SPRN = {{'SP1N','SP2N'}};
TLSE = {{'TL1E','TL2E','TL3E','TL4E'}};
TORI = {{'TO1I','TO2I','TO3I','TO4I'}};
TRF2 = {{'TR12'}};
UNPG = {{'UN1G','UN2G','UN3G','UN4G','UN5G'}};
VEN1 = {{'VE11'}};
VFCH = {{'VF1H','VF2H','VF3H'}};
WTZR = {{'WT1R','WT2R'}};
ZADA = {{'ZA1A'}};
ZIM2 = {{'ZI12','ZI22','ZI32'}};
ZIMM = {{'ZI1M','ZI2M'}};
ZOUF = {{'ZO1F'}};

SplittedSites = struct('AJAC',AJAC,'AUTN',AUTN,'AXPV',AXPV,'BRST',BRST, ...
           'BSCN',BSCN,'BZRG',BZRG,'CHIZ',CHIZ,'COMO',COMO,'DRES',DRES, ...
           'EGLT',EGLT,'ENTZ',ENTZ,'GARI',GARI,'GENO',GENO,'GRAS',GRAS, ...
           'GRAZ',GRAZ,'GSR1',GSR1,'IENG',IENG,'IGMI',IGMI,'KARL',KARL,'LINZ',LINZ, ...
           'LROC',LROC,'MARS',MARS,'MEDI',MEDI,'MLVL',MLVL,'MOPS',MOPS, ...
           'PFA2',PFA2,'PORE',PORE,'POTS',POTS,'PRAT',PRAT,'PUYV',PUYV, ...
           'ROVE',ROVE,'SBG2',SBG2,'SJDV',SJDV,'SPRN',SPRN,'TLSE',TLSE, ...
           'TORI',TORI,'TRF2',TRF2,'UNPG',UNPG,'VEN1',VEN1,'VFCH',VFCH, ...
           'WTZR',WTZR,'ZADA',ZADA,'ZIM2',ZIM2,'ZIMM',ZIMM,'ZOUF',ZOUF);
       
SplittedSites_list = {'AJAC','AUTN','AXPV','BRST','BSCN','BZRG','CHIZ','COMO','DRES','EGLT','ENTZ','GARI','GENO','GRAS','GRAZ','GSR1','IENG','IGMI','KARL','LINZ','LROC','MARS','MEDI','MLVL','MOPS','PFA2','PORE','POTS','PRAT','PUYV','ROVE','SBG2','SJDV','SPRN','TLSE','TORI','TRF2','UNPG','VEN1','VFCH','WTZR','ZADA','ZIM2','ZIMM','ZOUF'};

clear AJAC AUTN AXPV BRST BSCN BZRG CHIZ COMO DRES EGLT ENTZ GARI GENO GRAS GRAZ GSR1 IENG IGMI KARL LINZ LROC MARS MEDI MLVL MOPS PFA2 PORE POTS PRAT PUYV ROVE SBG2 SJDV SPRN TLSE TORI TRF2 UNPG VEN1 VFCH WTZR ZADA ZIM2 ZIMM ZOUF

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
            VEL_unified(row,:) = mean(VEL(iSet,:),1);
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
    VEL_ave(i,:) = mean(VEL_unified(iSet,:),1);
end

end