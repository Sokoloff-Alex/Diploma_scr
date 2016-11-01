% Data Processed

close all
clear all
clc

%%

fileNames = sort(listPSOL(:,1));
timeMartix_P = nan(10,366);
matrixWeeks_P = nan(7, 1773-1251);

for i = 1:length(fileNames)

    filename = cell2mat(fileNames(i));
    yy=str2num(filename(4:5));
    doy=str2num(filename(6:8));
    timeMartix_P(yy-3, doy) = 1;
    [status, stdout] = system(['gps_date -yd 20', filename(4:5), ' ', filename(6:8),' -o %W %w']);
    Week = str2num(stdout(1:4));
    dow = str2num(stdout(6));
    matrixWeeks_P(dow+1, Week-1250) = 1;

    
end
disp('Done')


%%

fileNames = sort(listSSOL(:,1));
timeMartix_S = nan(10,366);
matrixWeeks_S = nan(7, 1773-1251);
for i = 1:length(fileNames)
    filename = cell2mat(fileNames(i));
    yy=str2num(filename(4:5));
    doy=str2num(filename(6:8));
    timeMartix_S(yy-3, doy) = 1;
    [status, stdout] = system(['gps_date -yd 20', filename(4:5), ' ', filename(6:8),' -o %W %w']);
    Week = str2num(stdout(1:4));
    dow = str2num(stdout(6));
    matrixWeeks_S(dow+1, Week-1250) = 1;
    
end
disp('Done')


%% display missing
disp('list of missing sessions : ')
disp('yyyy ssss')

for iyy= 1:10
    for idd = 1:366
        if isnan(timeMartix_S(iyy,idd)) 
            
            if idd ~= 366 && ((iyy+2003) ~= 2005 || (iyy+2003) ~= 2006 || (iyy+2003) ~= 2007 || (iyy+2003) ~= 2009 || (iyy+2003) ~= 2010 || (iyy+2003) ~= 2011 || (iyy+2003) ~= 2013 || (iyy+2003) ~= 2014)
                disp([num2str(iyy+2003), ' ', num2str(idd),'0']);
            end
        end
    end
end


%%
close all
clc
figure(1)
subplot(4,1,1)
hold on
imagesc(timeMartix_P)
xlabel('doy')
ylabel('yy')
title('$P, green - available')
ax = gca;
set(ax,'YTick',[1:1:10])
set(ax,'YTickLabel',{[2004:1:2014]})
xlim([1 366])
ylim([1 10])


subplot(4,1,3)
hold on
imagesc(timeMartix_S)
xlabel('doy')
ylabel('yy')
title(' $SAVEDISK ,green - available')
ax = gca;
set(ax,'YTick',[1:1:10])
set(ax,'YTickLabel',{[2004:1:2014]})
xlim([1 366])
ylim([1 10])

subplot(4,1,2)
hold on
imagesc(matrixWeeks_P)
xlabel('week')
ylabel('dow')
PercentProcessed_P = nansum(nansum(matrixWeeks_P))/3653*100
title(['$P ,green - available : ', num2str(PercentProcessed_P), ' %'])
ax = gca;
set(ax,'YTick',[1:7])
set(ax,'YTickLabel',{[0:6]})
set(ax,'XTick',[1:52:522])
set(ax,'XTickLabel',{[1251:52:1773]})
xlim([1 1773-1251+1])
ylim([1 7])

subplot(4,1,4)
hold on
PercentProcessed_S = nansum(nansum(matrixWeeks_S))/3653*100
imagesc(matrixWeeks_S)
xlabel('week')
ylabel('dow')
title(['$SAVEDISK ,green - available : ', num2str(PercentProcessed_S), ' %'])
ax = gca;
set(ax,'YTick',[1:7])
set(ax,'YTickLabel',{[0:6]})
set(ax,'XTick',[1:52:522])
set(ax,'XTickLabel',{[1251:52:1773]})
% set(ax,'XTickLabel',{[2004:1:2014]})
set(ax,'XTickLabel',{[1251:52:1773]})
xlim([1 1773-1251+1])
ylim([1 7])

%%

get_Repeatability('Results/FCSIGSB.OUT')

%%
clc

listRINEXtimesort = load('Results/list_RINEX_time_sort');

%%
RINEX_Datapool = zeros(13,366);

filenames=listRINEXtimesort;

for i = 1:length(filenames)
    filename = cell2mat(filenames(i));
    doy=str2num(filename(5:7));
    yy=str2num(filename(10:11));
    RINEX_Datapool(yy-3, doy) = RINEX_Datapool(yy-3, doy) + 1;
end
disp('Done')

%% Asign RINEX to Networks

ALPEN      = {'AGNE','BOSC','BREI','CARZ','DEVE','ELMO','FAHR','FDOS','FERH','FERR','HELM','HGRA','HRIE','MAVE','MBEL','MITT','MOCA','OATO','PARO','POGG','PORA','SOND','WART'};
AUSTRIA    = {'GMND','GRAZ','HFL2','HFLK','HKBL','KOE2','KOET','KRBG','KTZ2','KTZB','PAT2','PATK','PFA2','PFAN','RIED','ROHR','SBG2','SBGZ','TRF2','TRFB','VLCH','VLKM','WELS','WIEN'};
FREDNET    = {'ACOM','AFAL','CANV','CODR','FUSE','JOAN','MDEA','MPRA','NOVE','PAZO','SUSE','TRIE','UDI1','UDIN','VARM','ZOUF'};
RENAG      = {'AGDE','AIGL','AJAC','ALPE','AUBU','AXPV','BANN','BART','BAUB','BUAN','CHAR','CHIZ','CHTL','CLFD','COMO','ENTZ','EOST','ERCK','ESAB','FJCP','GINA','GRAS','GROI','HEAU','JANU','JOUX','LACA','LROC','LUCE','MAKS','MANS','MARS','MICH','MODA','MOGN','MTPL','NICA','NICE','PALI','PARD','PQRL','PUYA','PUYV','RIXH','ROSD','SJDV','SLVT','STJ9','STMR','TLSE','TORI','VALD','WLBH'};
ORPHEON    = {'ALLE','ANNO','ARMI','AUMO','AVEN','BIWI','BLIX','BUIS','BURE','CAPA','CAZE','CHAM','CHMX','CHRN','CLAP','CURA','ESNO','FCLZ','FENO','FRAC','GUIL','JUVI','LAJA','LEBE','LFAZ','LUVI','MOLA','MOLV','MONT','MOUS','OGAG','PIGN','PLOE','POLI','PUEC','RABU','RSTL','SAPI','SAUV','SETE','SIMA','SOPH','SOUR','STEY','STGR','TENC','TRES','TROC','TROP','VAUC','VIGY'};
GREF       = {'DILL','ERLA','GOET','HOFJ','WT21'};
EUREF      = {'AUTN','AXPV','BOLG','BRST','BSCN','BZRG','CAME','CHIZ','COMO','DILL','DRES','EGLT','ENTZ','GARI','GENO','GSR1','HOFJ','IGMI','KARL','LINZ','MANS','MDOR','MLVL','MOPS','OBE2','OBE4','PFA2','PFAN','PORE','PRAT','PUYV','ROVE','SBG2','SBGZ','SJDV','SPRN','TORI','TRF2','TRFB','UNPG','VEN1','VENE','VFCH','ZADA','ZOUF'};
IGS        = {'AJAC','BRST','BZRG','FFMJ','GENO','GRAS','GRAZ','HFL2','HFLK','HUEG','IENG','LROC','MARS','MEDI','OBE2','OBE4','PADO','POTS','TLSE','VENE','WTZR','ZIM2','ZIMM'};

%

for i = 1:length(ALPEN);    Table1(i,1:2) = {ALPEN{i}  ,'ALPEN'};   end;
for i = 1:length(AUSTRIA);  Table2(i,1:2) = {AUSTRIA{i},'AUSTRIA'}; end;
for i = 1:length(FREDNET);  Table3(i,1:2) = {FREDNET{i},'FREDNET'}; end;
for i = 1:length(RENAG)     Table4(i,1:2) = {RENAG{i}  ,'RENAG'};   end;
for i = 1:length(ORPHEON)   Table5(i,1:2) = {ORPHEON{i},'ORPHEON'}; end;
for i = 1:length(GREF)      Table6(i,1:2) = {GREF{i},   'GREF'};    end;
for i = 1:length(EUREF)     Table7(i,1:2) = {EUREF{i},  'EUREF'};   end;
for i = 1:length(IGS)       Table8(i,1:2) = {IGS{i},    'IGS'};     end;

%

Table = [Table1; Table2; Table3; Table4; Table5; Table6; Table7; Table8]

Table_sorted = sortrows(Table)

%%
clc
tic

Matrix = zeros(13,366,9);
for i = 1:length(listRINEXtimesort)
   file = listRINEXtimesort(i);
   file = file{1};
   SITE = file(1:4);
   SITE = upper(SITE);
   yy =  str2num(file(10:11))-3;
   ddd = str2num(file(5:7));

   for j = 1:212;
        if strcmp(SITE,Table_sorted(j,1))
            Network = Table_sorted(j,2);   
            switch Network{:}
                case 'IGS';
                    net = 8;
                case 'EUREF'; 
                    net = 7;
                case 'ALPEN';
                    net = 1;
                case 'AUSTRIA';
                    net = 2;
                case 'FREDNET';
                    net = 3;
                case 'RENAG'; 
                    net = 4;
                case 'ORPHEON';
                    net = 5;
                case 'GREF';  
                    net = 6;
                otherwise 
                    disp(['Site ', mat2str(SITE) , ' in not assgned to the network'])
                    net = 9;
            end                      
            Matrix(yy,ddd,net) = Matrix(yy,ddd,net)+1;
            break;      
        end
   end
end
t1 = toc
disp('Done')
%%

ALPEN_rinex   = [Matrix(1,:,1)';Matrix(2,1:365,1)';Matrix(3,1:365,1)';Matrix(4,1:365,1)';Matrix(5,:,1)';Matrix(6,1:365,1)';Matrix(7,1:365,1)';Matrix(8,1:365,1)';Matrix(9,:,1)';Matrix(10,1:365,1)';Matrix(11,1:365,1)';Matrix(12,1:365,1)';Matrix(13,:,1)'];
AUSTRIA_rinex = [Matrix(1,:,2)';Matrix(2,1:365,2)';Matrix(3,1:365,2)';Matrix(4,1:365,2)';Matrix(5,:,2)';Matrix(6,1:365,2)';Matrix(7,1:365,2)';Matrix(8,1:365,2)';Matrix(9,:,2)';Matrix(10,1:365,2)';Matrix(11,1:365,2)';Matrix(12,1:365,2)';Matrix(13,:,2)'];
FREDNET_rinex = [Matrix(1,:,3)';Matrix(2,1:365,3)';Matrix(3,1:365,3)';Matrix(4,1:365,3)';Matrix(5,:,3)';Matrix(6,1:365,3)';Matrix(7,1:365,1)';Matrix(8,1:365,3)';Matrix(9,:,3)';Matrix(10,1:365,3)';Matrix(11,1:365,3)';Matrix(12,1:365,3)';Matrix(13,:,3)'];
RENAG_rinex   = [Matrix(1,:,4)';Matrix(2,1:365,4)';Matrix(3,1:365,4)';Matrix(4,1:365,4)';Matrix(5,:,4)';Matrix(6,1:365,4)';Matrix(7,1:365,4)';Matrix(8,1:365,4)';Matrix(9,:,4)';Matrix(10,1:365,4)';Matrix(11,1:365,4)';Matrix(12,1:365,4)';Matrix(13,:,4)'];
ORPHEON_rinex = [Matrix(1,:,5)';Matrix(2,1:365,5)';Matrix(3,1:365,5)';Matrix(4,1:365,5)';Matrix(5,:,5)';Matrix(6,1:365,5)';Matrix(7,1:365,5)';Matrix(8,1:365,5)';Matrix(9,:,5)';Matrix(10,1:365,5)';Matrix(11,1:365,5)';Matrix(12,1:365,5)';Matrix(13,:,5)'];
GREF_rinex    = [Matrix(1,:,6)';Matrix(2,1:365,6)';Matrix(3,1:365,6)';Matrix(4,1:365,6)';Matrix(5,:,6)';Matrix(6,1:365,6)';Matrix(7,1:365,6)';Matrix(8,1:365,6)';Matrix(9,:,6)';Matrix(10,1:365,6)';Matrix(11,1:365,6)';Matrix(12,1:365,6)';Matrix(13,:,6)'];
EUREF_rinex   = [Matrix(1,:,7)';Matrix(2,1:365,7)';Matrix(3,1:365,7)';Matrix(4,1:365,7)';Matrix(5,:,7)';Matrix(6,1:365,7)';Matrix(7,1:365,7)';Matrix(8,1:365,7)';Matrix(9,:,7)';Matrix(10,1:365,7)';Matrix(11,1:365,7)';Matrix(12,1:365,7)';Matrix(13,:,7)'];
IGS_rinex     = [Matrix(1,:,8)';Matrix(2,1:365,8)';Matrix(3,1:365,8)';Matrix(4,1:365,8)';Matrix(5,:,8)';Matrix(6,1:365,8)';Matrix(7,1:365,8)';Matrix(8,1:365,8)';Matrix(9,:,8)';Matrix(10,1:365,8)';Matrix(11,1:365,8)';Matrix(12,1:365,8)';Matrix(13,:,8)'];

%%
% # of sites estimated
clc
Sites_CRD_estimatedW = zeros(10,366);

filenames=listCRDestimatedW;

for i = 1:length(filenames)
    filename = cell2mat(filenames(i));
    doy=str2num(filename(15:17));
    yy=str2num(filename(13:14));
    Number=str2num(filename(24:end));
    Sites_CRD_estimatedW(yy-03, doy) = Number;
end
disp('Done')

Sites_CRD_estimatedA = zeros(10,366);

filenames=listCRDestimatedA;

for i = 1:length(filenames)
    filename = cell2mat(filenames(i));
    doy=str2num(filename(15:17));
    yy=str2num(filename(13:14));
    Number=str2num(filename(24:end));
    Sites_CRD_estimatedA(yy-03, doy) = Number;
end
disp('Done')

RINEX_computedW = [Sites_CRD_estimatedW(1,:),Sites_CRD_estimatedW(2,1:365),Sites_CRD_estimatedW(3,1:365),Sites_CRD_estimatedW(4,1:365),Sites_CRD_estimatedW(5,:),Sites_CRD_estimatedW(6,1:365),Sites_CRD_estimatedW(7,1:365),Sites_CRD_estimatedW(8,1:365),Sites_CRD_estimatedW(9,:),Sites_CRD_estimatedW(10,1:365)];
RINEX_computedA = [Sites_CRD_estimatedA(1,:),Sites_CRD_estimatedA(2,1:365),Sites_CRD_estimatedA(3,1:365),Sites_CRD_estimatedA(4,1:365),Sites_CRD_estimatedA(5,:),Sites_CRD_estimatedA(6,1:365),Sites_CRD_estimatedA(7,1:365),Sites_CRD_estimatedA(8,1:365),Sites_CRD_estimatedA(9,:),Sites_CRD_estimatedA(10,1:365)];



%%

RINEX_quantity = [RINEX_Datapool(1,:),RINEX_Datapool(2,1:365),RINEX_Datapool(3,1:365),RINEX_Datapool(4,1:365),RINEX_Datapool(5,:),RINEX_Datapool(6,1:365),RINEX_Datapool(7,1:365),RINEX_Datapool(8,1:365),RINEX_Datapool(9,:),RINEX_Datapool(10,1:365),RINEX_Datapool(11,1:365),RINEX_Datapool(12,1:365),RINEX_Datapool(13,:)];

clr = colormap(colorcube(10));

close all
figure(3)
hold on
grid on
startDate = datenum('01-01-2004');
endDate   = datenum('12-31-2016');
xData = linspace(startDate,endDate,365*13+4);
% plot(xData,RINEX_computedW,'r')
% plot(xData,RINEX_computedA,'g')
plot(xData,ALPEN_rinex,'.-','Color',clr(2,:))
plot(xData,AUSTRIA_rinex,'.-','Color',clr(3,:))
plot(xData,FREDNET_rinex,'.-','Color',clr(4,:))
plot(xData,RENAG_rinex,'.-','Color',clr(5,:))
plot(xData,ORPHEON_rinex,'.-','Color',clr(6,:))
plot(xData,GREF_rinex,'.-','Color',clr(7,:))
plot(xData,EUREF_rinex,'.-','Color',clr(8,:))
plot(xData,IGS_rinex,'.-','Color',clr(9,:))
% plot(xData,RINEX_quantity,'.-k')

ylabel(' # of RINEX files in $D/RINEX/')
datetick('x','yyyy')
legend('ALPEN','AUSTRIA','FREDNET','RENAG','ORPHEON','GREF','EUREF','IGS','Location','NorthEast' )


%%

close all
figure(3)
hold on
grid on
startDate = datenum('01-01-2004');
endDate   = datenum('12-31-2016');
xData = linspace(startDate,endDate,365*13+4);
plot(xData,FREDNET_rinex,'.-')
ylabel(' # of RINEX files in $D/RINEX/')
datetick('x','yyyy')

%%

[datestr(xData( (RENAG_rinex < 20 & RENAG_rinex > 0) ), 'yyyy-mm-dd '), num2str(RENAG_rinex (RENAG_rinex < 20 & RENAG_rinex > 0))]


%%

close all
figure(4)
hold on
grid on
startDate = datenum('01-01-2004');
endDate   = datenum('12-31-2013');
xData = linspace(startDate,endDate,365*10+3);
plot(xData,RINEX_computedW,'.-g')
plot(xData,RINEX_computedA,'.-m')
plot(xData,RINEX_quantity,'.-b')
plot(xData,RINEX_quantity-(RINEX_computedW + RINEX_computedA),'.-r')

ylabel(' # of RINEX files in $D/RINEX/')
datetick('x','yyyy')
legend('flag W','flag A','Total','Not processed','Location','NorthWest' )

%%
figure(4)
pcolor(RINEX_Datapool)



