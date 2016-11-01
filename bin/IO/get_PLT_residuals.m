function [meanResN, meanResE, meanResU, meanResH, meanRes3D, Sites, Domes] = get_PLT_residuals(filename, Station_list)
% parse PLT file and compute residuals RMS for list of station
% 
% Alexandr Sokolov, KEG
% 12.10.2016

%%

filename_tmp = 'tmp.PLT';

disp('SITE DOMES_nmr         RMS_N,mm     RMS_E,mm     RMS_U,mm     RMS_H,mm     RMS_T,mm')
len = length(Station_list);

Sites = repmat(char(0),len,4);
Domes = repmat(char(0),len,9);
meanResN  = zeros(len,1);
meanResE  = zeros(len,1);
meanResU  = zeros(len,1);
meanResH  = zeros(len,1);
meanRes3D = zeros(len,1);

for i = 1:len
    SITE = Station_list{i};
%     disp(['grep ''', SITE, ''' ', filename , ' > ' , filename_tmp]);
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
        resN    = str2double(lineN(25:35));
        resE    = str2double(lineE(25:35));
        resU    = str2double(lineU(25:35));   
        res_N_stack( counter,1) = resN;
        res_E_stack( counter,1) = resE;
        res_U_stack( counter,1) = resU;
        res_H_stack( counter,1) = sqrt(resN^2 + resE^2);
        res_3D_stack(counter,1) = norm([resN, resE, resU]);
        status = feof(fileID);
    end    
    fclose(fileID);
    Sites(i,:)  = lineN(1:4);
    Domes(i,:)  = lineN(6:14);

    meanResN(i)  = rms(res_N_stack);
    meanResE(i)  = rms(res_E_stack);
    meanResU(i)  = rms(res_U_stack);
    meanResH(i)  = rms(res_H_stack);
    meanRes3D(i) = rms(res_3D_stack);
    
    fprintf('%3d %4s %9s  %12.3f %12.3f %12.3f %12.3f %12.3f   \n',i, Sites(i,:),Domes(i,:), meanResN(i)*1000, meanResE(i)*1000, meanResU(i)*1000, meanResH(i)*1000, meanRes3D(i)*1000)
    else
        disp([SITE,' is not found in ', filename]) 
    end
    delete(filename_tmp);
end