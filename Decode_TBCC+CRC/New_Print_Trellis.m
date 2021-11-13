function [Printed]=New_Print_Trellis(Trellis_matrix, Trellis_state, Trellis_input_matrix, filename, k, m, list_size, g_TBCC, g_CRC)

% Note that list size is valid only for P-LVA, while use the value 2 for S-LVA

Printed = 0;
num_states=size(Trellis_matrix,1);
num_input_bits=log2(max(Trellis_input_matrix,[],'all')+1);
num_output_bits=size(Trellis_matrix,2) / size(Trellis_state,2);
num_next_states=size(Trellis_state,2);

fileID=fopen(filename,'w');
fprintf(fileID, '#include <stdint.h>\n\n#define k %d\n#define n %d\n#define m %d\n#define NumSections %d\n#define NumOutputs %d\n#define NumStates %d\n#define NumInputs %d\n#define NumNextStates %d\n#define r %d\n#define L %d\n#define Memory %d\n\n', k, num_output_bits*k, m, k+m, num_output_bits, num_states, num_input_bits, num_next_states, m, list_size, log2(num_states));

fprintf(fileID, 'typedef struct connections_node_t_s {\n// static info of the node\nint16_t next_states[%d];\n// array with indexes of the next states of this node (e.g. this node is connected to node 0 and 5, then next_states[0]=0 and next_states[1]=5)\nint16_t previous_states[%d];\n// array with indexes of the previous states of this node (e.g. this node is connected to node 0 and 5, then previous_states[0]=0 and previous_states[1]=5)\ndouble outputs[%d][%d];\n// matrix where every row contains the output symbols (NOT bits) from this node to the next state stored in that row index (e.g. output symbols when go to state 5 are (+1,-1) then outputs[1][0]=+1 and outputs[1][1] = -1)\nuint8_t inputs[%d][%d];\n// matrix where every row contains the input bits to generate the corresponding output symbols from this node to the next state stored in that row index (e.g. output symbols when go to state 5 are (+1,-1) and the input bit to go from this node to the next node is 1, then inputs[1][0]=1)\nuint8_t previous_states_inputs[%d][%d];\n// matrix where every row contains the input bits to generate the corresponding output symbols from the previous node in that row index and this node (e.g. output symbols when go to state 5 are (+1,-1) and the input bit to go from this node to the next node is 1, then inputs[1][0]=1)\n} connections_node_t;\n\n',num_next_states,num_next_states,num_next_states,num_output_bits,num_next_states,num_input_bits,num_next_states,num_input_bits);

fprintf(fileID,'typedef struct decoder_node_t_s {\ndouble edge_costs[%d];\n// array with the euclidean cost to go from this state to the corresponding next_states in that index\n// info of the local list of paths stored in this node. The first path is the one at minimum distance, while the second entry is the path from the other incoming node.\n// paths info : (list index, starting node , previous node , list index of the path at previous node , distance from the received vector)\nint16_t starting_node[%d];\n// array with indexes of the starting nodes of the stored paths\nint16_t previous_node[%d];\n// array with indexes of the previous nodes of the stored paths\nint16_t previous_index[%d];\n// array with list indexes of the stored paths at previous node list\ndouble path_metric[%d];\n// path metrics of the stored paths up to this node\n} decoder_node_t;\n\n',num_next_states,list_size,list_size,list_size,list_size);

if list_size >= 2
    fprintf(fileID,'typedef struct small_decoder_node_t_s {\ndouble edge_costs[%d];\n// array with the euclidean cost to go from this state to the corresponding next_states in that index\n// info of the local list of paths stored in this node. The first path is the one at minimum distance, while the second entry is the path from the other incoming node.\n// paths info : (list index, starting node , previous node , list index of the path at previous node , distance from the received vector)\nint16_t starting_node[%d];\n// array with indexes of the starting nodes of the stored paths\nint16_t previous_node[%d];\n// array with indexes of the previous nodes of the stored paths\nint16_t previous_index[%d];\n// array with list indexes of the stored paths at previous node list\ndouble path_metric[%d];\n// path metrics of the stored paths up to this node\n} small_decoder_node_t;\n\n',num_next_states,1,1,1,1);
end

fprintf(fileID, 'extern void initialize_trellis(connections_node_t (*Trellis_links)[NumStates]);');
fclose(fileID);

fileID = fopen(strcat(filename(1:end-1),'c'),'w');
fprintf(fileID, '#include "%s"\n\n',filename);
fprintf(fileID, 'extern uint8_t g_TBCC[NumOutputs][Memory+1]={');
for i_poly = 1:size(g_TBCC,1)
    fprintf(fileID, '{ %d', g_TBCC(i_poly,1));
    for i_bit = 2:size(g_TBCC,2)
    	fprintf(fileID, ', %d', g_TBCC(i_poly,i_bit));
    end
    if i_poly<size(g_TBCC,1)
	fprintf(fileID, '},');
    else
    fprintf(fileID, '}};\n');
    end
end
fprintf(fileID, 'extern uint8_t g_CRC[m+1]=');
fprintf(fileID, '{ %d', g_CRC(1));
for i_bit = 2:size(g_CRC,2)
fprintf(fileID, ', %d', g_CRC(i_bit));
end
fprintf(fileID, '};\n\n');

fprintf(fileID, 'void initialize_trellis(connections_node_t (*Trellis_links)[NumStates]){\n\n');
for i_state=1:num_states
    [predecessors,~] = find(Trellis_state==i_state);
    for i_next_state=1:num_next_states
    	bits=str2num(dec2bin(Trellis_input_matrix(i_state, Trellis_state(i_state,i_next_state)),num_input_bits).');
        bits_predecessor=str2num(dec2bin(Trellis_input_matrix(predecessors(i_next_state), i_state),num_input_bits).');
	fprintf(fileID,'    (*Trellis_links)[%d].next_states[%d]=%d;\n',i_state-1,i_next_state-1,Trellis_state(i_state,i_next_state)-1);
	fprintf(fileID,'    (*Trellis_links)[%d].previous_states[%d]=%d;\n',i_state-1,i_next_state-1,predecessors(i_next_state)-1);
	for i_out = 1:num_output_bits
		fprintf(fileID,'    (*Trellis_links)[%d].outputs[%d][%d]=%.1f;\n',i_state-1,i_next_state-1,i_out-1,Trellis_matrix(i_state,(i_next_state-1)*num_next_states+i_out));
	end
	for i_in = 1:num_input_bits
		fprintf(fileID,'    (*Trellis_links)[%d].inputs[%d][%d]=%d;\n',i_state-1,i_next_state-1,i_in-1,bits(i_in));
	end
	for i_in = 1:num_input_bits
		fprintf(fileID,'    (*Trellis_links)[%d].previous_states_inputs[%d][%d]=%d;\n',i_state-1,i_next_state-1,i_in-1,bits_predecessor(i_in));
	end
    end
end

fprintf(fileID, '\n}\n');
fclose(fileID);
Printed=1;






end
