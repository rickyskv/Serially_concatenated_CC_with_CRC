// Serial List (Wrap-Around) Viterbi Algorithm

// Note: this file contains a simple implementation of the S-LVA using the tree-trellis list Viterbi algorithm (see Roder and Hamzaoui, 2006).
// It is not multi-threading. The multiple processes are used to speed-up simulations, where each process is an indipendent decoder/simulator.
// A central processor receives from the other processes the number of codewords simulated before reaching a number of errors defined by the user.
// To set properly the parameters of this C code, use the associated MATLAB file.

// Possible improvements will be published in the future:
// - use 2 stacks to save memory in the decoder
// - quantize the decoder
// - use intel vector instructions to speed-up the forward passages when fixed-point arithmetic is used and also in backward
// - substitute the v = u G encoder with the convolutional encoder
// - substitute the syndrome check via matrix multiplication with the one via LSFR used by the BCH/CRC codes
// - improve central process communication abilities to avoid strugglers

// ##################################
// --- libraries

//libraries for vector operations (not used)
#include "emmintrin.h"
#include "immintrin.h"
#include "xmmintrin.h"
#include "tmmintrin.h"

//libraries for multi-processes in LINUX
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
#include <time.h>

//standard libraries
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
//#include <malloc.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "./trellis.h"

// --- definitions
#define NAME_SIZE 100
#define PI 3.1415926536
#define START "A"

// ############################################
// --- import global variables 
uint8_t g_CRC[m+1];
uint8_t g_TBCC[NumOutputs][Memory+1];

// ######################################
// --- structures

// structure of the stack entries used to store the paths in S-LVA
typedef struct stack_node_t_s {
  int16_t starting_node; // to check easily the tail-biting condition
  int16_t ending_node; // to check easily the tail-biting condition
  int16_t *section_split_indexes; // section indexes at which this path stops following the best path and it splits to the other branch
  int16_t *state_split_indexes; // state indexes where this path goes when it splits to the other branch
  int16_t num_splits; // number of splits
  double path_metric; // path metric of that path
  double partial_path_metric; // path metric from split to end 
} stack_node_t;

// heap data-structure used to keep the stack entries sorted
typedef struct heap_t_s{
    stack_node_t *nodes;
    int len;
    int size;
} heap_t;

// #############################################
// --- headers of the functions
// The functions are reported after the main

// convert dB in linear
double db2dec(double A){
    return pow(10.0,A/10.0);
}
// TBCC encoder
void TBCC_encode (uint8_t message[k+m], uint8_t (*codeword)[NumSections*NumOutputs]);
// CRC encoder
void CRC_encode (uint8_t message[k], uint8_t (*codeword)[k+m]);
// CRC check
int CRC_check (uint8_t message[k+m]);
// insert a new node in the heap datastructure
void heap_push (heap_t *h, stack_node_t new_node);
// extract the best node in the heap datastructure
stack_node_t heap_pop (heap_t *h);
// read parity check matrix
void WAVA_Update_Distances_and_Build_Tree_Trellis(connections_node_t Trellis_links[NumStates], decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], uint8_t num_WAVA);
// Serial List Viterbi Algorithm
void Serial_List_Viterbi_Decoder(decoder_node_t Trellis_nodes[NumSections+1][NumStates], int list_size, int *needed_list_size, uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK);
// backward function to find the codeword/message corresponding to a certain path stored in the stack
void back_Trellis_stack(decoder_node_t Trellis_nodes[NumSections+1][NumStates],stack_node_t stack_node, uint8_t (*codeword)[NumSections*NumInputs]);
// Additive-White Gaussian Noise generator
double AWGN_generator();
// update the trellis at every new codeword decoding
void update_trellis (connections_node_t Trellis_links[NumStates], decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], double Y[NumSections*NumOutputs], int list_size);
// function to set communication pipes between processes
void set_pipe(int *** fd, int num_proc);//set number of pipes to use
// read a puncturing pattern
void read_puncturing(char *Puncturing_file, uint8_t (*puncturing_v)[NumSections*NumOutputs]);
// functions to free dynamical arrays
void free2D(int **Matrix, int num_rows);
void free2D_unsigned(unsigned int **Matrix, int num_rows);


// ###########################################################
// --- MAIN FUNCTION
int main(int argc, char *argv[])
{
    
    char *output_filename, *buf, *puncturing_vector_filename; // buf is used to communicate between processes

    int num_errors_per_process, count_errors, count_UE, count_NACK, tot_NACK, num_errors, i_check, i_codeword, i_bit, stop, count_punct, punct_every, num_punctures, num_processes, i_process, **f_pipe, **q_pipe, num_codewords_per_process, status,i_u_bit, total_bits, i_row, actual_list_size, total_list_size; // k is the number of info bits, n the number of codeword bits of the convolutional code without CRC, m is the degree of the CRC code, ML_message is the found maximum likelihood message (if anything is found), r is the number of rows of the parity check matrix, i_ are indexes used, stop is a flag, f_pipe and q_pipe are used by the processes to communicate, status is used to check the status of the processes
	
    uint8_t NACK, ML_message[NumSections * NumInputs], num_WAVA_iterations;
    int16_t list_size;
	
    uint8_t u[k], u_CRC[k+m], v[NumSections * NumOutputs], all_zero[NumSections * NumOutputs], puncturing_v[NumSections * NumOutputs]; // u contains the info bit sent (message), v is the codeword sent, G the generator matrix of the code, all_zero is a all_zero vector
    
    decoder_node_t Trellis_nodes[NumSections+1][NumStates]; // Trellis structure
    connections_node_t Trellis_links[NumStates]; // Trellis structure
    
    double Y[NumSections * NumOutputs], start_EbN0, end_EbN0, step_EbN0, CER, i_EbN0, std_noise, temp_double_value, average_list_size; // Y is the modulated codeword, all EbN0 terms are in dB, CER is the Codeword Error Rate, std_noise the standard deviation of the noise, while temp_double_value is used as temporary variable

    FILE *file_out; // file pointer
    pid_t pid; // process ID for the multi-process
    
    // CHECK THAT THERE ARE ALL ARGUMENTS AND SAVE THEM, OTHERWISE RETURN AN ERROR
    if (argc != 10){
        printf("ERROR! The program arguments shall be:\n./PLVA output_filename num_WAVA_iterations list_size num_errors puncting_vector_filename num_processes start_EbN0 end_EbN0 step_EbN0\n");
        exit(1);
    }
    else {
        output_filename=argv[1];
        num_WAVA_iterations=atoi(argv[2]); // if 1, pure Viterbi
        list_size=atoi(argv[3]);
        num_errors=atoi(argv[4]); // number of errors to be counted at each E_b/N_0 simulated point
        puncturing_vector_filename = argv[5];
    	num_processes=atoi(argv[6]);
        start_EbN0=atof(argv[7]);
        end_EbN0=atof(argv[8]);
        step_EbN0=atof(argv[9]);
    }
    
    

    initialize_trellis(&Trellis_links);

     //printf("Out of the function values :%d, %d, %.1f\n",Trellis_nodes[2][0].next_states[0],Trellis_nodes[2][0].next_states[1],Trellis_nodes[0][0].outputs[1][1]);

    

 	// read puncturing vector and count punctures
 	read_puncturing(puncturing_vector_filename, &puncturing_v);

	num_punctures = 0;
	for (i_bit = 0; i_bit < NumSections * NumOutputs; i_bit++) num_punctures += (1-puncturing_v[i_bit]);
    
    // initialize pipes for communications among processes
    set_pipe(&f_pipe, num_processes);//on f_pipe children wait
    set_pipe(&q_pipe, num_processes);//on q_pipe children send results to the parent

    // set number of errors to be computed by each process
    num_errors_per_process = (int)ceil((double)num_errors/(double)num_processes); 

    // initialize the all-zero vector and the vector with the modulated symbols
    for(i_bit = 0; i_bit < NumSections * NumOutputs; i_bit++) all_zero[i_bit] = 0;
    for(i_bit=0;i_bit <NumSections*NumOutputs;i_bit++) Y[i_bit]= 0.0;


	total_bits = NumSections*NumOutputs;
    
    // generate the children processes
    for(i_process=0;i_process<num_processes;i_process++){
        //SET communications
        if(pipe(f_pipe[i_process])!=0 || pipe(q_pipe[i_process])!=0){
            fprintf(stderr,"ERROR with pipe\n");
            exit(99);
        }
        // fork a new process
        pid=fork();
        if(pid<0){//ERROR with fork()
            fprintf(stderr,"ERROR while creating new process!\n");
            exit(99);
        }
        else{
            if(pid==0){//CHILD		
                close(f_pipe[i_process][1]);//READ from f_pipe
                close(q_pipe[i_process][0]);//WRITE in q_pipe

                while(read(f_pipe[i_process][0],&buf[i_process],1)>0){}//wait until all children are created
                close(f_pipe[i_process][0]);
                srand((i_process+1)*time(0)); // set different random seed for each child process
                
                for (i_row = 0; i_row < k; i_row++) u[i_row]=0;
                for (i_bit = 0; i_bit < total_bits; i_bit++) v[i_bit]=0;
                
                // for loop over the different values of E_b / N_0
                for (i_EbN0 = start_EbN0; i_EbN0 <= end_EbN0; i_EbN0 += step_EbN0){
                    // initialize variables and compute noise standard deviation
                    CER = 0.0;
                    total_list_size=0;
                    std_noise=sqrt(1.0/(2.0*((double)k/(((double)(k+m)*(double)NumOutputs/(double)NumInputs)-(double)num_punctures)*db2dec(i_EbN0))));
                    
                    count_errors = 0;
                    count_UE = 0;
                    count_NACK = 0;
                    i_codeword = 0;

			for (i_row = 0; i_row < k; i_row++) u[i_row]=0;
					
                    // while loop of encoding -> channel -> decoding until num_errors_per_process errors are counted by the process
            		while (count_errors < num_errors_per_process){

                        // generate random message u and encode it via G to obtain the transmitted codeword v
			for (i_bit=0;i_bit<k;i_bit++) u[i_bit]=((uint8_t)rand()) % 2;	
			CRC_encode (u, &u_CRC);
			TBCC_encode (u_CRC, &v);
						
                        // generate transmitted modulated codeword via BSPK mapping (0 -> +1, 1 -> -1), puncture it (0->0,1->0) and add noise
                        count_punct=0;
                        for(i_bit=0;i_bit<total_bits;i_bit++){
                            if (puncturing_v[i_bit]==0){
                            	Y[i_bit]= (double)0.0;
                            }
                            else {
                            	(v[i_bit] == 0) ? (Y[i_bit]= 1.0 + (std_noise * AWGN_generator())) : (Y[i_bit]= -1.0 + (std_noise * AWGN_generator()));
                            }
                        }

                        // update the weigths of the edges of the trellis with the distance from the received codeword
                        update_trellis (Trellis_links, &Trellis_nodes, Y, 2) ;
                        

                        // if more than one WAVA iteration, run WAVA with list 1 to update the initial distances of the nodes in the initial trellis section. At the last WAVA iteration build the tree-trellis
                        WAVA_Update_Distances_and_Build_Tree_Trellis(Trellis_links, &Trellis_nodes, num_WAVA_iterations);

                        // run SLVA
                        Serial_List_Viterbi_Decoder(&Trellis_nodes, list_size, &actual_list_size, &ML_message, &NACK);
                        
                        // update expected list size measurement
                        total_list_size+= actual_list_size;
                        
                        // check if the found message is correct, otherwise increment the counter of the errors

                        if (NACK==0){
                                stop=0;
		                for(i_bit=0;i_bit<NumSections*NumInputs-m && stop==0; i_bit++){
					if (ML_message[i_bit] != u[i_bit]){
						stop=1;
						count_errors++;
						count_UE++;
					}
		                }
		        }
		        else{
		            count_NACK++;
		            count_errors++;
		        }
                        
                        // free previously allocated memory
						
                        // increment number of transmitted codewords
            		i_codeword++;

                    } // end while loop over the counter number of erroneus received codewords

                    // send to the parent process the number of transmitted codewords to reach num_errors_per_process errors
                    if(write(q_pipe[i_process][1],&i_codeword,sizeof(int))!=sizeof(int)){
                        fprintf(stderr,"ERROR while writing on pipe\n");
                        exit(99);
					}
		    // send to the parent process the total list size used
		    if(write(q_pipe[i_process][1],&total_list_size,sizeof(int))!=sizeof(int)){
                        fprintf(stderr,"ERROR while writing on pipe\n");
                        exit(99);
					}
		    // send to the parent process the number of NACKs
		    if(write(q_pipe[i_process][1],&count_NACK,sizeof(int))!=sizeof(int)){
                        fprintf(stderr,"ERROR while writing on pipe\n");
                        exit(99);
					}
                } // end for loops over E_b/N_0
                
                // close pipes to communicate with the parent and kill the child process
                close(q_pipe[i_process][1]);
                
                // free allocated memory
                free2D(f_pipe,num_processes);
                free2D(q_pipe,num_processes);

                
		        exit(EXIT_SUCCESS);
            }
	
            else{ //PARENT process
                // close not used pipes to communicate
                close(f_pipe[i_process][0]);//WRITE on f_pipe
                close(q_pipe[i_process][1]);//READ on q_pipe

                // send message to children processes to start their simulation
                if(write(f_pipe[i_process][1],START,1)!=1){
                    fprintf(stderr,"ERROR while writing on pipe\n");
                    exit(99);
                }
            }
        }
    } // close for loop over the processes generation
    
    // only parent process access this part of the main function
    
    for(i_process=0; i_process<num_processes; i_process++){//close f_pipe to let the different processes to start parallel computations
        close(f_pipe[i_process][1]);
    }

    // open the file where to write the codeword error rate results
    file_out = fopen(output_filename, "w");
    if (!file_out) {
        printf("ERROR! Problems in writing in the file : <%s>\n", output_filename);
        exit(12);
    }

    // for loop over E_b/N_0 to collect the results computed by children processes
    for (i_EbN0 = start_EbN0; i_EbN0 <= end_EbN0; i_EbN0 += step_EbN0){
        CER = 0.0;
        tot_NACK = 0;
        average_list_size = 0.0;
        // collect the result of each process and sum them
        for(i_process=0; i_process<num_processes; i_process++){
            read(q_pipe[i_process][0],&num_codewords_per_process,sizeof(int));
            read(q_pipe[i_process][0],&total_list_size,sizeof(int));
            read(q_pipe[i_process][0],&count_NACK,sizeof(int));
            CER += (double)num_codewords_per_process; // count total simulated codewords
            average_list_size += (double)total_list_size;
            tot_NACK +=  count_NACK;
        }
        average_list_size /= CER;

        CER = (double)(num_processes*num_errors_per_process) / CER;

        // print result for each E_b / N_0 in the file and in the standard output
        //fprintf(file_out,"%f\n",CER);
        fprintf(file_out,"Eb/N0 = %.2f , CER = %.9f, average L = %.2f, NACK/(NACK+UE) %d/%d\n",i_EbN0, CER, average_list_size, tot_NACK,num_processes*num_errors_per_process);
        printf("Eb/N0 = %.2f , CER = %.9f, average L = %.2f, NACK/(NACK+UE) %d/%d\n",i_EbN0, CER, average_list_size, tot_NACK,num_processes*num_errors_per_process);
    }

    // when processed all E_b / N_0 close the output file and wait that each child process is terminated correctly
    fclose(file_out);
    for(i_process=0; i_process<num_processes; i_process++){//read on q_pipe
        while(read(q_pipe[i_process][0],&temp_double_value,sizeof(double))>0){}
        close(q_pipe[i_process][0]);
        wait(&status);
    }

    // free allocated memory

    free2D(f_pipe,num_processes);
    free2D(q_pipe,num_processes);
    
    // return that all is ok
    return 0;
}

// ##################################
// --- functions description



// update the weights of the edges of the trellis structure with the one of the new received vector
void update_trellis (connections_node_t Trellis_links[NumStates], decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], double Y[NumSections*NumOutputs], int list_size){

    // define variables

    int16_t i_section, i_state, i, j;
    double temp_edge_cost, temp_value;
    
    // compute the new weights of each edge and store them in the structure
    for (i_section = 0; i_section < NumSections; i_section ++){
        for(i_state = 0; i_state < NumStates; i_state++){
            for (i=0 ; i< NumNextStates; i++){
                temp_edge_cost = 0;
                for (j=0; j< NumOutputs; j++){
                    temp_value=Y[i_section*NumOutputs+j] - Trellis_links[i_state].outputs[i][j];
                    temp_edge_cost += temp_value*temp_value;
                }
                (*Trellis_nodes)[i_section][i_state].edge_costs[i] = temp_edge_cost;
            }
    	}
    }
            // update list infos
        
    for(i_state = 0; i_state < NumStates; i_state++){
	(*Trellis_nodes)[0][i_state].starting_node[0] = i_state;
	(*Trellis_nodes)[0][i_state].path_metric[0] = 0.0;
	(*Trellis_nodes)[0][i_state].starting_node[1] = -1;
	(*Trellis_nodes)[0][i_state].path_metric[1] = -1.0;
    }
}

/*TBCC encoder using the convolutional encoder g_TBCC = (g_0, g_1, ..., g_v), v is the memory*/
void TBCC_encode (uint8_t message[k+m], uint8_t (*codeword)[NumSections*NumOutputs]){
	uint8_t state[Memory];
	int i_state, i_out, i_bit;

	i_bit=k+m-1;
	// initialize finite state machine
	for (i_state = 0; i_state < Memory; i_state++) state[i_state] = message[i_bit-i_state];
	
	// initialize codeword to all 0
	for (i_bit = 0; i_bit < NumSections*NumOutputs; i_bit++) (*codeword)[i_bit]=0;
	
	// encode using convolutional encoder
	for (i_bit=0; i_bit<k+m;i_bit++){		
		for (i_out=0; i_out<NumOutputs; i_out++){
			// encode
			for(i_state=0; i_state<Memory; i_state++) (*codeword)[i_bit*NumOutputs+i_out]^=(g_TBCC[i_out][i_state+1] & state[i_state]);//state[Memory-i_state-1]);
			(*codeword)[i_bit*NumOutputs+i_out]^=(g_TBCC[i_out][0] & message[i_bit]);//(g_TBCC[i_out][Memory] & message[i_bit]);
		}
		// update state
		for(i_state=Memory-1; i_state>0; i_state--) state[i_state] = state[i_state-1];
		state[0] = message[i_bit];
	}
	
}

/*CRC systematic encoder using the LSFR encoder g_CRC = (g_0, g_1, ..., g_v)*/
void CRC_encode (uint8_t message[k], uint8_t (*codeword)[k+m]){
	uint8_t state[m], new_bit;
	int i_state, i_bit;
	
	// initialize finite state machine
	
	for (i_state = 0; i_state < m; i_state++) state[i_state]=0;
	
	// systematic -> copy first k bits, while updating the state
	for (i_bit = 0; i_bit < k; i_bit++){
		(*codeword)[i_bit]=message[i_bit]^0;
		new_bit = message[i_bit] ^ state[m-1];
		
		for (i_state=m-1;i_state>0;i_state--) state[i_state] = (new_bit & g_CRC[i_state])^state[i_state-1];
		state[0] = new_bit;
	}
	// copy last bits
	i_bit = k;
	for (i_state = m-1; i_state >=0; i_state--) (*codeword)[i_bit+(m-1-i_state)] = state[i_state];
}

/*CRC systematic encoder using the LSFR encoder g_CRC = (g_0, g_1, ..., g_v)*/
int CRC_check (uint8_t message[k+m]){
	uint8_t codeword[k+m], partial_message[k];
	
	int i_bit;
	
	for (i_bit = 0; i_bit < k; i_bit++) partial_message[i_bit]=message[i_bit];
	
	//CRC_encode (g_CRC, partial_message, &codeword);
	CRC_encode (partial_message, &codeword);	
	
	for (i_bit = k; i_bit < k+m; i_bit++){
		if (codeword[i_bit]!=message[i_bit]) return 0;
	}
	
	return 1;
}


/*
wrap-around Viterbi decoder for the first (num_WAVA - 1) iterations of the WAVA algorithm, to set the initial distances of the nodes in first section
plus this function builds the tree-trellis structure: each node has a list of 2 made of the best incoming path in the first entry, plus path from the other incoming edge in the second entry
*/
void WAVA_Update_Distances_and_Build_Tree_Trellis(connections_node_t Trellis_links[NumStates], decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], uint8_t num_WAVA){
    // declare variables
    int16_t j_state, index[NumStates][2], bit[NumStates][2], i_WAVA, i_value[NumStates], min_index[NumStates];
    double metric[NumStates][2], other_path_metric;
    register int16_t i_state;
    int16_t i_section;
	
    // set to 0 the path metric at each node in the last trellis section
    for (i_state = 0; i_state < NumStates ; i_state++) (*Trellis_nodes)[NumSections][i_state].path_metric[0]=0.0;
    
    // for loop over the WAVA iterations
    for (i_WAVA=0 ; i_WAVA < num_WAVA ; i_WAVA++){
        // update every state in the first trellis section at every new WAVA iteration, with the values of the same nodes in the final trellis section 
        for (i_state = 0; i_state < NumStates ; i_state++){
            (*Trellis_nodes)[0][i_state].starting_node[0] = i_state;
            (*Trellis_nodes)[0][i_state].previous_node[0] = -1;
            (*Trellis_nodes)[0][i_state].previous_index[0] = -1;
            (*Trellis_nodes)[0][i_state].path_metric[0] = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
        }
        
	// run Viterbi and store both best path and other incoming path to each node at each trellis section 
        for (i_section=1 ; i_section <= NumSections ; i_section++){
				
		for (i_state = 0; i_state < NumStates ; i_state++){
			// extract data from the 2 predecessors
			index[i_state][0] = Trellis_links[i_state].previous_states[0];
			index[i_state][1] = Trellis_links[i_state].previous_states[1];
			bit[i_state][0] = Trellis_links[i_state].previous_states_inputs[0][0];
			bit[i_state][1] = Trellis_links[i_state].previous_states_inputs[1][0];
			metric[i_state][0] = (*Trellis_nodes)[i_section-1][index[i_state][0]].path_metric[0] + (*Trellis_nodes)[i_section-1][index[i_state][0]].edge_costs[bit[i_state][0]];
			metric[i_state][1] = (*Trellis_nodes)[i_section-1][index[i_state][1]].path_metric[0] + (*Trellis_nodes)[i_section-1][index[i_state][1]].edge_costs[bit[i_state][1]];
			
			// check which has the best path
			(metric[i_state][0] <= metric[i_state][1]) ? (min_index[i_state] = 0) : (min_index[i_state] = 1);
			
			// store best
			(*Trellis_nodes)[i_section][i_state].starting_node[0] = (*Trellis_nodes)[i_section-1][index[i_state][min_index[i_state]]].starting_node[0];
			(*Trellis_nodes)[i_section][i_state].previous_node[0] = index[i_state][min_index[i_state]];
			(*Trellis_nodes)[i_section][i_state].previous_index[0] = bit[i_state][min_index[i_state]];
			(*Trellis_nodes)[i_section][i_state].path_metric[0] = metric[i_state][min_index[i_state]];
			
			// store other
			(*Trellis_nodes)[i_section][i_state].starting_node[1] = (*Trellis_nodes)[i_section-1][index[i_state][1-min_index[i_state]]].starting_node[0];
			(*Trellis_nodes)[i_section][i_state].previous_node[1] = index[i_state][1-min_index[i_state]];
			(*Trellis_nodes)[i_section][i_state].previous_index[1] = bit[i_state][1-min_index[i_state]];
			(*Trellis_nodes)[i_section][i_state].path_metric[1] = metric[i_state][1-min_index[i_state]];
			
		} // close state for loop
        } // close section for loop
        
		
    } // close num_WAVA for loop

}

/*
Serial List Viterbi Algorithm run at the last WAVA iteration
*/
void Serial_List_Viterbi_Decoder(decoder_node_t Trellis_nodes[NumSections+1][NumStates], int list_size, int *needed_list_size, uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK){
    // define variables
    int16_t i_bit, i_section, i_state, i_list_state, stop = 0, best_codeword_index, best_state, counter, j_state, index, last_section, found_a_codeword;
    int i_list, actual_list_size, next_list_size, i_stack;
    uint8_t temp_codeword[NumSections*NumInputs];

	
	heap_t *heap_stack=(heap_t *)calloc(1, sizeof (heap_t));

	stack_node_t stack_node, best_node;
    	int stack_size=0,temp_stack1_size=0, temp_stack2_size=0, list_index=0, new_split, count_list ;
    	int16_t starting_section, starting_state;
	double worst_path_metric, temp_path_metric, partial_path_metric;

    found_a_codeword = 0;
    actual_list_size = 1;
    next_list_size = 1;
		
	// store in the heap datastructure all the paths found at the nodes in the last trellis section
    // remember that each node in the stack represents a path and it has the following data:
    // (starting node, ending node, path metric, number of splits from the best path, section indexes of the splits, nodes of the splits, partial path metric from the last split to the ending node)
	stack_size=0;
	worst_path_metric = -1.0;
    for (i_state = 0; i_state < NumStates ; i_state++){
		stack_node.starting_node = Trellis_nodes[NumSections][i_state].starting_node[0];
		stack_node.ending_node = i_state;
		stack_node.num_splits = 0;
		stack_node.path_metric = Trellis_nodes[NumSections][i_state].path_metric[0];
		stack_node.partial_path_metric=0.0;
		
		if (worst_path_metric < stack_node.path_metric) worst_path_metric = stack_node.path_metric;
		
        // insert in the heap datastructure
		heap_push(heap_stack,stack_node);
		stack_size++;
	}
	
			
	// check if best path is tail-biting and it matches the CRC condition. If this is the case, stop
	// extract node from the stack
	best_node = heap_pop(heap_stack);
	stack_size--;
    // check tail-biting condition
	if (best_node.starting_node == best_node.ending_node){
        // extract message of the corresponding path
		back_Trellis_stack(Trellis_nodes,best_node,&temp_codeword);
		//check CRC condition via syndrome check
		if (CRC_check(temp_codeword)==1){
            	// save as best message
			found_a_codeword = 1;
			for (i_bit=0;i_bit<NumSections*NumInputs;i_bit++) (*best_codeword)[i_bit]=temp_codeword[i_bit];
		}
	}
	
	// counter of the list size
	count_list = 1;
	
	// if first entry is not a codeword of the TBCC+CRC code:
    // 1- go backward in that path and add the local best paths to that path in the stack as new nodes. use the heap to keep the stack sorted
    // 2- remove that node from the stack
    // 3- extract new best node in the stack and check tail-biting condition and CRC condition
    //    a- if both conditions are met stop
    //    b- otherwise repeat from 1- to 3- until you find a path respecting both conditions or you have checked L paths, with L the maximum allowed list size
	while (found_a_codeword == 0 && count_list < list_size && stack_size>0){
		// backward + add to stack the local best paths to that path
		
        // check if best path has splits or not
		if (best_node.num_splits > 0) { // if there are splits, you start going backward from the last split
            // set starting section
			starting_section = best_node.section_split_indexes[best_node.num_splits-1]-1;
            // set starting node
			starting_state = best_node.state_split_indexes[best_node.num_splits-1];
			i_state = Trellis_nodes[starting_section+1][starting_state].previous_index[1];
			starting_state = Trellis_nodes[starting_section+1][starting_state].previous_node[1];
            // initialize the partial path metric (euclidean distance of the path from the last split to the ending node w.r.t. the received vector)
			partial_path_metric = best_node.partial_path_metric + Trellis_nodes[starting_section][starting_state].edge_costs[i_state];
            i_state = starting_state;
		}
		else{ // if there are no splits you start going backward from the last trellis section
            // set starting section
			starting_section = NumSections;
            // set starting node
			starting_state = best_node.ending_node;
			i_state = starting_state;
            // initialize the partial path metric (euclidean distance of the path from the last split to the ending node w.r.t. the received vector). Note that in thiscase is 0.0
			partial_path_metric = best_node.partial_path_metric;
		}
		
		// set number of splits
		new_split = best_node.num_splits;
		
        // go backward from last split (if any)
		for(i_section = starting_section; i_section > 0; i_section--){
			
            // compute path metric of the local best path
			temp_path_metric = Trellis_nodes[i_section][i_state].path_metric[1] + partial_path_metric;
			
            // insert in the list if you have not reached L entries or, in case the list is full, insert it if it is better than the last entry
            // to be modified later with a max-heap structure keeping the best L (max list value) path metrics and using the L-th path metric as threshold for insertion
			if (temp_path_metric < worst_path_metric || stack_size < list_size){
                // copy info of the path in a new stack entry
				stack_node.starting_node = Trellis_nodes[i_section][i_state].starting_node[1];
				stack_node.ending_node = best_node.ending_node;
				stack_node.num_splits = best_node.num_splits+1;
				
                // copy the splits of the path
				stack_node.section_split_indexes=(int*)malloc(stack_node.num_splits*sizeof(int));
				stack_node.state_split_indexes=(int*)malloc(stack_node.num_splits*sizeof(int));
					
				if (stack_node.section_split_indexes == NULL || stack_node.state_split_indexes==NULL){
					printf("error while allocating for splits of the stack node!\n");
					exit(118);
				}
				
				if (best_node.num_splits > 0){
					memcpy(stack_node.section_split_indexes,best_node.section_split_indexes,best_node.num_splits*sizeof(int));
					memcpy(stack_node.state_split_indexes,best_node.state_split_indexes,best_node.num_splits*sizeof(int));
				}
				
				stack_node.section_split_indexes[new_split] = i_section;
				stack_node.state_split_indexes[new_split] = i_state;
				
				stack_node.path_metric = temp_path_metric;
				stack_node.partial_path_metric = partial_path_metric;
				
				if (worst_path_metric < stack_node.path_metric) worst_path_metric = stack_node.path_metric;
                
                // insert the new stack entry in the stack in a sorted position
				heap_push(heap_stack,stack_node);

                // increase stack size
				stack_size++;
                
                // deallocate not used memory
				free(stack_node.section_split_indexes);
				free(stack_node.state_split_indexes);
			}
			
			// update partial path metric adding the edge distance to go to the node in the previous section of the actual best path. Remember we are going backward
			if (i_section > 0) partial_path_metric += Trellis_nodes[i_section-1][Trellis_nodes[i_section][i_state].previous_node[0]].edge_costs[Trellis_nodes[i_section][i_state].previous_index[0]];
			
            // find node index of the actual best path at the previous section
			i_state = Trellis_nodes[i_section][i_state].previous_node[0];
			
			
		}// end extraction of local best paths to the actual best path
		
		// deallocate not used memory
		if(best_node.num_splits>0){
			free(best_node.section_split_indexes);
			free(best_node.state_split_indexes);
		}
		
        // pick next best path and check tail-biting and CRC conditions, if ok STOP
		best_node = heap_pop(heap_stack);
		stack_size--;
		// check tail-biting
		if (best_node.starting_node == best_node.ending_node){
			back_Trellis_stack(Trellis_nodes,best_node,&temp_codeword);
			// check CRC syndrome
			if (CRC_check(temp_codeword)==1){
				found_a_codeword = 1;
   				for (i_bit=0;i_bit<NumSections*NumInputs;i_bit++) (*best_codeword)[i_bit]=temp_codeword[i_bit];
				
			}
		}
        
        // increase list counter. Rememebr we can check up to L paths
		count_list++;
	} // end while loop to find a tail-biting path which meets the CRC condition among the best L trellis paths
	
	// deallocate memory of the stack
	for (i_stack = 0; i_stack < stack_size;i_stack++){
		stack_node=heap_pop(heap_stack);
		if(stack_node.num_splits>0){
			free(stack_node.state_split_indexes);
			free(stack_node.section_split_indexes);
		}
	}
	free(heap_stack->nodes);
	free(heap_stack);
	
	// if you have not found a path which respects the tail-biting and CRC constraints, then the algorithm outputs an erasure (a message of all -1 in this implementation)
    (found_a_codeword == 0) ? (*NACK = 1) : (*NACK=0);
    //if (found_a_codeword == 0){
	//for (i_section = 0; i_section < NumSections*NumInputs ; i_section++) (*best_codeword)[i_section] = -1;
    //}
    
    *needed_list_size = count_list;
    
}

// go bacward in the trellis and extract the input message associated to a path stored in the stack
void back_Trellis_stack(decoder_node_t Trellis_nodes[NumSections+1][NumStates],stack_node_t stack_node, uint8_t (*codeword)[NumInputs*NumSections]){

    int16_t i_section, i_input, i_next_state, counter_splits, initial_state, previous_state, actual_state, previous_list_index, i_bit, stop=0;


    initial_state = stack_node.ending_node;
    actual_state = initial_state;

    i_bit=NumSections - 1;
	
	
	if (stack_node.num_splits == 0) {
		for (i_section = NumSections; i_section > 0 ; i_section--){
			(*codeword)[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[0];
			i_bit --;
			actual_state = Trellis_nodes[i_section][actual_state].previous_node[0];
		}
	}
	else {
		counter_splits = 0;
		for (i_section = NumSections; i_section > 0 ; i_section--){
			if (i_section == stack_node.section_split_indexes[counter_splits]){
				(*codeword)[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[1];
				i_bit --;
				actual_state = Trellis_nodes[i_section][actual_state].previous_node[1];
				if (counter_splits  < stack_node.num_splits) counter_splits++;
			}
			else{
				(*codeword)[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[0];
				i_bit --;
				actual_state = Trellis_nodes[i_section][actual_state].previous_node[0];
			}
		}
		
	}
	
}



// read a puncturing pattern
void read_puncturing(char *Puncturing_file, uint8_t (*puncturing_v)[NumSections*NumOutputs]){
    /*
    The file containing the puncturing vector must be of the form :
    
        p_0 p_1 p_2 ... p_{n-1}
    
    with p_i = 0 if punctured bit, or 1 if NOT punctured bit
    */
    FILE *file_P;
    int i_bit;
    unsigned int N=NumSections*NumOutputs, temp; 
    
    file_P = fopen(Puncturing_file, "r");
    if (file_P) {  
        for(i_bit=0;i_bit<N; i_bit++){
	    fscanf(file_P, "%d",&temp);
	    (*puncturing_v)[i_bit]=(uint8_t) temp;
        }
        fclose(file_P);
    }
    else{
        printf("ERROR! Problems in reading the file : <%s>\n", Puncturing_file);
        exit(15);
    }
    
}

// from https://www.embeddedrelated.com/showcode/311.php
double AWGN_generator()
{/* Generates additive white Gaussian Noise samples with zero mean and a standard deviation of 1. */
  double temp1;
  double temp2;
  double result;
  int p;
  p = 1;
  while( p > 0 )
  {
	temp2 = ( rand() / ( (double)RAND_MAX ) ); /*  rand() function generates an integer between 0 and  RAND_MAX, which is defined in stdlib.h. */
    if ( temp2 == 0 )
    {// temp2 is >= (RAND_MAX / 2)
      p = 1;
    }// end if
    else
    {// temp2 is < (RAND_MAX / 2)
       p = -1;
    }// end else
  }// end while()
  temp1 = cos( ( 2.0 * (double)PI ) * rand() / ( (double)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;
  return result;	// return the generated random sample to the caller
}// end AWGN_generator()

// function to set the pipes to communicate among processes
void set_pipe(int *** fd, int num_proc){
	int i, **f_pipe;

	f_pipe=(int **)malloc(num_proc*sizeof(int *));

	if(f_pipe==NULL){
		fprintf(stderr,"ERROR while allocating memory\n");
		exit(99);
	}

	for(i=0; i<num_proc; i++){
		f_pipe[i]=(int *)malloc(2*sizeof(int));
		if(f_pipe[i]==NULL){
			fprintf(stderr,"ERROR while allocating memory\n");
			exit(99);
		}	
	}

	*fd=f_pipe;
}

// functions to manage the heap datastructure
// insert a new node
void heap_push (heap_t *h, stack_node_t new_node) {
	
    if (h->len + 1 >= h->size) {
        h->size = h->size ? h->size * 2 : 4;
        h->nodes = (stack_node_t *)realloc(h->nodes, h->size * sizeof (stack_node_t));
    }
    int i = h->len + 1;
    int j = i / 2;
    while (i > 1 && h->nodes[j].path_metric > new_node.path_metric) {
        h->nodes[i] = h->nodes[j];
        i = j;
        j = j / 2;
    }
    h->nodes[i].starting_node = new_node.starting_node;
    h->nodes[i].ending_node = new_node.ending_node;

    if(new_node.num_splits > 0){
		h->nodes[i].section_split_indexes = (int *)malloc(new_node.num_splits*sizeof(int));
		h->nodes[i].state_split_indexes = (int *)malloc(new_node.num_splits*sizeof(int));
		
		memcpy(h->nodes[i].section_split_indexes, new_node.section_split_indexes, new_node.num_splits*sizeof(int));
		memcpy(h->nodes[i].state_split_indexes, new_node.state_split_indexes, new_node.num_splits*sizeof(int));
	}

    h->nodes[i].num_splits = new_node.num_splits;
    h->nodes[i].path_metric = new_node.path_metric;
    h->nodes[i].partial_path_metric = new_node.partial_path_metric;

    h->len++;
}
// extract the best node 
stack_node_t heap_pop (heap_t *h) {
    int i, j, kk;
    
	
    stack_node_t best_node = h->nodes[1];
	
    h->nodes[1] = h->nodes[h->len];
 
    h->len--;
 
    i = 1;
    while (i!=h->len+1) {
        kk = h->len+1;
        j = 2 * i;
        if (j <= h->len && h->nodes[j].path_metric < h->nodes[kk].path_metric) {
            kk = j;
        }
        if (j + 1 <= h->len && h->nodes[j + 1].path_metric < h->nodes[kk].path_metric) {
            kk = j + 1;
        }
        h->nodes[i] = h->nodes[kk];
        i = kk;
    }
	
    return best_node;
}

// functions to deallocate memory
void free2D(int **Matrix, int num_rows){
    int i_row;
    
    for(i_row=0; i_row<num_rows; i_row++) free(Matrix[i_row]);
    free(Matrix);
}

void free2D_unsigned(unsigned int **Matrix, int num_rows){
    int i_row;
    
    for(i_row=0; i_row<num_rows; i_row++) free(Matrix[i_row]);
    free(Matrix);
}
