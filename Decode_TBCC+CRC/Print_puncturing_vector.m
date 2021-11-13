function [Printed]=Print_puncturing_vector(punct_v, filename)

Printed = 0;
n= numel(punct_v);

fileID=fopen(filename,'w');
for i_bit=1:n
    fprintf(fileID,'%d ',punct_v(i_bit));
end
fclose(fileID);

Printed=1;
