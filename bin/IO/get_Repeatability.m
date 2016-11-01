function [Records] = get_Repeatability(filename);
%% parse  *.OUT file, extract repeatability tables and plot them
    close all
    clear all
    clc

  filename = 'Results/FCSIGSB.OUT';

%% extract part of file 
    disp(['Parsing ', filename, ' ...']);
    [status, StartLine]= system(['grep --line-number "Comparison of individual solutions:"  ',filename,' | cut -f1 -d:']);
    [status, EndLine]=   system(['grep --line-number "Variance-covariance scaling factors:" ',filename,' | cut -f1 -d:']);
    Len = str2num(EndLine) - str2num(StartLine);

    [status] = system(['tail -n +', num2str(str2num(StartLine)), ' ', filename,'| head -',num2str(Len),' > Results/Repeatability.txt' ]);

    [status, NumberOfStations] =   system(['grep -c " U " Results/Repeatability.txt']);
    NumberOfStations = str2num(NumberOfStations);
    [status, NumberOf_NEQ_Files] = system(['grep -c "/SOL/" ', filename]);
    NumberOfWeeks = str2num(NumberOf_NEQ_Files)/7;

    filename = 'Results/Repeatability.txt';
  %%
    fileID = fopen(filename, 'r');

    disp(fgets(fileID));
    disp(fgets(fileID));
    disp(fgets(fileID));

    StationNamesStack = [];
    ResidualsMatrix = nan(NumberOfStations,3,NumberOfWeeks+1,7);
%%
for iStation = 1:NumberOfStations
    for iComponent = 1:3
        tline_head = fgets(fileID); 
        Station = tline_head([2:5, 7:15]);
        if length(tline_head)>28; 
            iMax = (length(tline)-30)/8;
            for iDow = 1:iMax
                Residual = str2num(tline_head((30+(iDow*8)-7):(30+(iDow*8)))); %#ok<ST2NM>
                if ~isempty(Residual)
                    ResidualsMatrix(iStation,iComponent, 1, iDow) = Residual;
                end
            end
        end
        for iRow = 2:NumberOfWeeks
            tline = fgets(fileID);
            if length(tline)>1; 
               iMax = (length(tline)-30)/8;
               for iDow = 1:iMax
                    Residual = str2double(tline((30+(iDow*8)-7):(30+(iDow*8))));
                    if ~isempty(Residual)
                        ResidualsMatrix(iStation,iComponent, iRow, iDow) = Residual;
                    end
               end
            end
        end
        
    end
%     disp(tline_head)
    tline = fgets(fileID);
    StationNamesStack = [StationNamesStack; Station];
    StationRepeatability = struct('N', squeeze(ResidualsMatrix(iStation,1,:,:)), 'E', squeeze(ResidualsMatrix(iStation,2,:,:)), 'U', squeeze(ResidualsMatrix(iStation,3,:,:)));
    
    if iStation == 1
        Records = struct(Station, StationRepeatability);
    else
        Records = setfield(Records, StationNamesStack(end,:), StationRepeatability);
    end
    
end

fclose(fileID);
clear ResidualsMatrix StationRepeatability tline tline_head A 
clear i iRow iStation iDow iMax iComponent counter Station ans Residual 

disp('Done')

%%
fields = fieldnames(Records);
B = [];
for iOld = 1:numel(fields)
    A = [Records.(fields{iOld}).E;  Records.(fields{iOld}).N; Records.(fields{iOld}).U];
    B = [B,  A];
end

%%
clc

fields = fieldnames(Records);
C = [];
StationNumber = numel(fields);

iOld = 1;
site = fields{1};
site_old = site(1:4);
domes_old = site(5:end);
buffer = [Records.(fields{1}).E;  Records.(fields{1}).N; Records.(fields{1}).U];
C = [];

for iOld = 2:NumberOfStations
    site = fields{iOld};
    site_new= site(1:4);
    domes_new = site(5:end);
    new = [Records.(fields{iOld}).E; Records.(fields{iOld}).N; Records.(fields{iOld}).U];
    
    if strcmp(site_old, site_new) %%|| strcmp(domes_old, domes_new)
        buff(1,:,:) = buffer;
        buff(2,:,:) = new;
        buffer = squeeze(nansum(buff,1));
    else
        C = [C buffer];
        buffer = new;
    end  
    site_old = site_new;
    domes_old = domes_new;
end

C = [C new];

%%
figure(2)
image(B)
colorbar

figure(3)
image(C)
colorbar

figure(4)
imagesc(B)
colorbar
caxis([-10 10])

figure(5)
imagesc(C)
colorbar
caxis([-10 10])

%%

figure(4)
hold on
grid on
plot(nanmean(B, 2), 'b')
plot(nanmean(C, 2), 'r')

%%
[status] =   system(['rm Results/Repeatability.txt']);
disp('Done')

end
