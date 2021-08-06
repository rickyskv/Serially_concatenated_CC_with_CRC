function [Printed]=Print_G(G, filename)

Printed = 0;
r=size(G,1);
n=size(G,2);

fileID=fopen(filename,'w');
for i_row=1:r
    for i_column=1:n
        fprintf(fileID,'%d ',G(i_row,i_column));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);

Printed=1;