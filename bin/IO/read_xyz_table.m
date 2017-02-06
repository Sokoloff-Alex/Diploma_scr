function [out_array] = read_xyz_table(filename)
% function to read table from file with 3 columns

[status, msgout]= system(['cat ', filename, ' | wc --lines ']);
NumOfLines = str2double(msgout);

out_array = zeros(NumOfLines,4);

FileID = fopen(filename);
for i = 1:NumOfLines
    line = fgetl(FileID);
    out_array(i,1) = str2double(line(1:2));  % yy
    out_array(i,2) = str2double(line(3:5));  % ddd
    out_array(i,3) = str2double(line(7:10)); % Num of baselines
    out_array(i,4) = str2double(line(12:16));% Resolved (in %)
end

fclose(FileID);
end