function [Printed]=Print_CRC(H, filename)

Printed = 0;
r=size(H,1);
n=size(H,2);

fileID=fopen(filename,'w');
fprintf(fileID,'%d\n',r);

for i_row=1:r
    for i_column=1:n
        fprintf(fileID,'%d ',H(i_row,i_column));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);

Printed=1;