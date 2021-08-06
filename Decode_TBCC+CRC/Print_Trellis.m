function [Printed]=Print_Trellis(Trellis_matrix, Trellis_state, Trellis_input_matrix, filename)

Printed = 0;
num_states=size(Trellis_matrix,1);
num_input_bits=log2(max(Trellis_input_matrix,[],'all')+1);
num_output_bits=size(Trellis_matrix,2) / size(Trellis_state,2);
num_next_states=size(Trellis_state,2);

fileID=fopen(filename,'w');
fprintf(fileID,'%d %d %d\n',num_states,...
    num_input_bits, ...
    num_next_states);

for i_state=1:num_states
    for i_next_state=1:num_next_states
        fprintf(fileID,'%d ',Trellis_state(i_state,i_next_state)-1);
    end
    for i_output_bit=1:num_output_bits*num_next_states
        fprintf(fileID,'%.1f ',Trellis_matrix(i_state,i_output_bit));
    end
    for i_next_state=1:num_next_states
        bits=str2num(dec2bin(Trellis_input_matrix(i_state, Trellis_state(i_state,i_next_state)),num_input_bits).');
        for i_input_bit=1:numel(bits)
            fprintf(fileID,'%d ',bits(i_input_bit));
        end
    end
    fprintf(fileID,'\n');
end
fclose(fileID);
Printed=1;