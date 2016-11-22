function [Stations, Radoms, Records] = readOUT(filename)
% Parse *.OUT file from Bernese   
% exctract only crd/vel information table
%
% Alexandr Sokolo, KEG
% 2016

tmpFile = 'tmp.txt'; 

if exist(filename, 'file') ~= 2
    disp(['ERROR, File not found: ',filename])
else
    fullpath = which(filename);


% Extract part of file 


if exist(fullpath, 'file')
    % File exists.  Do stuff....
    %filename = 'Results/FCSIGSB.OUT';
    disp(['Parsing ', fullpath, ' ...']);
    [status, StartLine]= system(['grep --line-number "Station name          Typ   A priori value  Estimated value    Correction     RMS error      3-D ellipsoid        2-D ellipse"  ',filename,' | cut -f1 -d:']);
    [status, EndLine]  = system(['grep --line-number "Helmert Transformation Parameters With Respect to Combined Solution:" ',fullpath,' | cut -f1 -d:']);
    Len = str2num(EndLine) - str2num(StartLine)-2;
    [status] = system(['tail -n +', num2str(str2num(StartLine)), ' ', fullpath,'| head -',num2str(Len),' > ' , tmpFile ]);
else
  % File does not exist.
  warningMessage = sprintf('Warning: file does not exist:\n%s', fullpath);
%   uiwait(msgbox(warningMessage));
end


   

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*1s%4s%1s%9s%8s%1s%19s%17s%14s%14s%12s%7s%12s%s%[^\n\r]';

%% Open the text file.

fileID = fopen(tmpFile,'r');
disp(fgetl(fileID));
disp(fgetl(fileID));
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[4,5,6,7,8,9,10,11]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

% %% Extract StationNames
% clc
% 
% NameColumn = dataArray{1};
% Radoms = dataArray{3};
% 
% NumberOfStations = size(NameColumn,1)/12;  % 4 entry type x 3 component
% for i = 1:NumberOfStations
%     Station(i) = NameColumn(i*12-11,:);
%     Radom(i) = Radoms(i*12-11,:);
% end


%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [4,5,6,7,8,9,10,11]);
rawCellColumns = raw(:, [1,2,3]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Clear temporary variables
clearvars filename formatSpec fileID dataArray ans col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;

[status, cmdout] = system(['rm ', tmpFile]);

%% Rearrange table

% NumberOfStations = size(raw,1)/12;

    Stations = raw(1:12:end,1);
%     Stations = cell2mat(Stations);
%     Stations = cellstr(Stations);
    Radoms   = raw(1:12:end,3);
    Radoms = cell2mat(Radoms);
%     NumberOfStations = length(Stations);
    
%     Coordinates
    XYZ_Apr = [cell2mat(raw(1:12:end,6)), cell2mat(raw(2:12:end,6)) ,cell2mat(raw(3:12:end,6))]; 
    XYZ_Est = [cell2mat(raw(1:12:end,7)), cell2mat(raw(2:12:end,7)) ,cell2mat(raw(3:12:end,7))]; 
    XYZ_Cor = [cell2mat(raw(1:12:end,8)), cell2mat(raw(2:12:end,8)) ,cell2mat(raw(3:12:end,8))]; 
    XYZ_rms = [cell2mat(raw(1:12:end,9)), cell2mat(raw(2:12:end,9)) ,cell2mat(raw(3:12:end,9))]; 
    
    ENU_Apr = [cell2mat(raw(6:12:end,6)), cell2mat(raw(5:12:end,6)) ,cell2mat(raw(4:12:end,6))]; 
    ENU_Est = [cell2mat(raw(6:12:end,7)), cell2mat(raw(5:12:end,7)) ,cell2mat(raw(4:12:end,7))]; 
    ENU_Cor = [cell2mat(raw(6:12:end,8)), cell2mat(raw(5:12:end,8)) ,cell2mat(raw(4:12:end,8))]; 
    ENU_rms = [cell2mat(raw(6:12:end,9)), cell2mat(raw(5:12:end,9)) ,cell2mat(raw(4:12:end,9))]; 
        
    ENU_3D_ellipsoid_sigmas = [cell2mat(raw(6:12:end,10)), cell2mat(raw(5:12:end,10)) ,cell2mat(raw(4:12:end,10))];
    ENU_3D_ellipsoid_angles = [cell2mat(raw(6:12:end,11)), cell2mat(raw(5:12:end,11)) ,cell2mat(raw(4:12:end,11))];
    EN_2D_ellipse_sigmas = [str2num(cell2mat(raw(6:12:end,12))), str2num(cell2mat(raw(5:12:end,12)))];
    EN_2D_ellipse_azimuth = str2num(cell2mat(raw(5:12:end,13)));
    
    %     Velocities
    V_XYZ_Apr = [cell2mat(raw(7:12:end,6)), cell2mat(raw(8:12:end,6)) ,cell2mat(raw(9:12:end,6))]; 
    V_XYZ_Est = [cell2mat(raw(7:12:end,7)), cell2mat(raw(8:12:end,7)) ,cell2mat(raw(9:12:end,7))]; 
    V_XYZ_Cor = [cell2mat(raw(7:12:end,8)), cell2mat(raw(8:12:end,8)) ,cell2mat(raw(9:12:end,8))]; 
    V_XYZ_rms = [cell2mat(raw(7:12:end,9)), cell2mat(raw(8:12:end,9)) ,cell2mat(raw(9:12:end,9))]; 
    
    V_ENU_Apr = [cell2mat(raw(12:12:end,6)), cell2mat(raw(11:12:end,6)) ,cell2mat(raw(10:12:end,6))]; 
    V_ENU_Est = [cell2mat(raw(12:12:end,7)), cell2mat(raw(11:12:end,7)) ,cell2mat(raw(10:12:end,7))]; 
    V_ENU_Cor = [cell2mat(raw(12:12:end,8)), cell2mat(raw(11:12:end,8)) ,cell2mat(raw(10:12:end,8))]; 
    V_ENU_rms = [cell2mat(raw(12:12:end,9)), cell2mat(raw(11:12:end,9)) ,cell2mat(raw(10:12:end,9))]; 
    V_ENU_3D_ellipsoid_sigmas = [cell2mat(raw(12:12:end,10)), cell2mat(raw(11:12:end,10)) ,cell2mat(raw(10:12:end,10))];
    V_ENU_3D_ellipsoid_angles = [cell2mat(raw(12:12:end,11)), cell2mat(raw(11:12:end,11)) ,cell2mat(raw(10:12:end,11))];
    V_EN_2D_ellipse_sigmas = [str2num(cell2mat(raw(12:12:end,12))), str2num(cell2mat(raw(11:12:end,12))) ]; 
    V_EN_2D_ellipse_azimuth = str2num(cell2mat(raw(11:12:end,13)));
    
    %%
    XYZ = struct('Apr', XYZ_Apr,'Est', XYZ_Est, 'Cor',XYZ_Cor,'rms', XYZ_rms);
    Ellipsoid = struct('Sigmas', ENU_3D_ellipsoid_sigmas, 'Angles', ENU_3D_ellipsoid_angles);
    Ellipse = struct('Sigmas', EN_2D_ellipse_sigmas, 'Angle', EN_2D_ellipse_azimuth);
    ENU = struct('Apr', ENU_Apr, 'Est', ENU_Est, 'Cor', ENU_Cor, 'rms', ENU_rms, 'Ellipsoid', Ellipsoid,'Ellipse', Ellipse);
    CRD = struct('XYZ', XYZ, 'ENU', ENU);
    
    V_XYZ = struct('Apr', V_XYZ_Apr, 'Est', V_XYZ_Est, 'Cor',V_XYZ_Cor, 'rms',V_XYZ_rms);
    V_Ellipsoid = struct('Sigmas', V_ENU_3D_ellipsoid_sigmas, 'Angles', V_ENU_3D_ellipsoid_angles);
    V_Ellipse = struct('Sigmas', V_EN_2D_ellipse_sigmas, 'Angle', V_EN_2D_ellipse_azimuth);
    V_ENU = struct('Apr', V_ENU_Apr, 'Est', V_ENU_Est, 'Cor',V_ENU_Cor, 'rms',V_ENU_rms, 'Ellipsoid',V_Ellipsoid, 'Ellipse',V_Ellipse);
    VEL = struct('XYZ',V_XYZ, 'ENU', V_ENU);  
    
    
%     Records = struct('CRD', CRD, 'VEL', VEL, 'NumberOfStations', NumberOfStations, 'Stations', Stations, 'Radoms', Radoms);
    Records = struct('CRD', CRD, 'VEL', VEL, 'Stations', cell2mat(Stations), 'Radoms', Radoms);
    
    disp('Done')
    
end
     
end
    
    
