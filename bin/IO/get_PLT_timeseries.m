function [res_N_ts, res_E_ts, res_U_ts, MJD, epochNmr] = get_PLT_timeseries(filename, SITE)
% parse PLT file for timeseries for a single SITE 
% 
% Alexandr Sokolov, KEG
% 12.10.2016

%%
filename_tmp = 'tmp.PLT';
[status_nix]= system(['grep ''', SITE, ''' ', filename , ' > ' , filename_tmp]);
if status_nix == 0
    fileID = fopen(filename_tmp,'r');
    status = feof(fileID);
    res_N  = nan(12.5*366,1);
    res_E  = nan(12.5*366,1);
    res_U  = nan(12.5*366,1);    
    res_H  = nan(12.5*366,1);
    res_3D = nan(12.5*366,1);
    counter = 0;
    while ~status
        counter = counter + 1;
        lineN = fgetl(fileID);
        lineE = fgetl(fileID);
        lineU = fgetl(fileID);
        res_N_ts( counter,1) = str2double(lineN(25:35));
        res_E_ts( counter,1) = str2double(lineE(25:35));
        res_U_ts( counter,1) = str2double(lineU(25:35));     
        epochNmr(counter,1)  = str2num(lineN(18:21));
        MJD(counter,1)       = str2num(lineN(36:48));
        status = feof(fileID);
    end
    fclose(fileID);
    delete(filename_tmp);
    Site     = lineN(1:4);
    DomesNmr = lineN(6:14);
else
    disp([SITE,' is not found in ', filename]) 
end
end