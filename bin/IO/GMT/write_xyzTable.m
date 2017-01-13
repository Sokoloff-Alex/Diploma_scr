function write_xyzTable(Records, filename, formatStr)
% write SIMPLE XYZ table for GMT 

fileID = fopen(filename, 'w');

disp([' ... wrinting xyz table into ', filename])

for i = 1:size(Records,1)    
    fprintf(fileID, formatStr, Records(i,:));   
end
disp('Done')
fclose(fileID);

end