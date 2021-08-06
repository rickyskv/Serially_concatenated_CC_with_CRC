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
#include <malloc.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

// --- definitions
#define NAME_SIZE 100
#define PI 3.1415926536
#define START "A"

// ######################################
// --- structures

// trellis node structure
typedef struct decoder_node_t_s {
  // static info of the node
  int *next_states; // array with indexes of the next states of this node (e.g. this node is connected to node 0 and 5, then next_states[0]=0 and next_states[1]=5)
  double **outputs; // matrix where every row contains the output symbols (NOT bits) from this node to the next state stored in that row index (e.g. output symbols when go to state 5 are (+1,-1) then outputs[1][0]=+1 and outputs[1][1] = -1)
  int **inputs; // matrix where every row contains the input bits to generate the corresponding output symbols from this node to the next state stored in that row index (e.g. output symbols when go to state 5 are (+1,-1) and the input bit to go from this node to the next node is 1, then inputs[1][0]=1)
  double *edge_costs; // array with the euclidean cost to go from this state to the corresponding next_states in that index
  
  // info of the local list of paths stored in this node. The first path is the one at minimum distance, while the second entry is the path from the other incoming node.
  // paths info : (list index, starting node , previous node , list index of the path at previous node , distance from the received vector) 
  int starting_node[2]; // array with indexes of the starting nodes of the stored paths
  int previous_node[2]; // array with indexes of the previous nodes of the stored paths
  int previous_index[2]; // array with list indexes of the stored paths at previous node list
  double path_metric[2]; // path metrics of the stored paths up to this node
} decoder_node_t;

// structure of the stack entries used to store the paths in S-LVA
typedef struct stack_node_t_s {
  int starting_node; // to check easily the tail-biting condition
  int ending_node; // to check easily the tail-biting condition
  int *section_split_indexes; // section indexes at which this path stops following the best path and it splits to the other branch
  int *state_split_indexes; // state indexes where this path goes when it splits to the other branch
  int num_splits; // number of splits
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
// insert a new node in the heap datastructure
void heap_push (heap_t *h, stack_node_t new_node);
// extract the best node in the heap datastructure
stack_node_t heap_pop (heap_t *h);
// read parity check matrix
void read_H (char *H_filename, int ***H, int n, int *r);
// build the trellis structure
decoder_node_t **initialize_trellis (char *trellis_filename, double *Y, int num_sections, int *num_states, int num_outputs, int list_size, int *num_inputs, int *num_next_states) ;
// run Viterbi to update the initial distance of each starting node of the trellis and build the tree trellis structure (each node stores the best path according to Viterbi and the path from the other incoming node) 
void WAVA_Update_Distances_and_Build_Tree_Trellis(decoder_node_t ***Trellis_nodes, int num_sections, int num_states, int num_next_states, int num_WAVA);
int *Serial_List_Viterbi_Decoder(decoder_node_t ***Trellis_nodes, int num_sections, int num_states, int num_next_states, int list_size, int num_inputs, int **H, int r);
// Parallel List Viterbi Algorithm
int *Parallel_List_Viterbi_Decoder(decoder_node_t ***Trellis_nodes, int num_sections, int num_states, int num_next_states, int list_size, int num_inputs, int **H, int r);
// backward function to find the codeword/message corresponding to a certain path stored in the stack
int *back_Trellis_stack(decoder_node_t **Trellis_nodes,stack_node_t stack_node,int num_sections,int num_next_states, int num_inputs);
// CRC syndrome check of the paths found
int check_syndrome(int *temp_codeword, int **H, int n, int r);
// Additive-White Gaussian Noise generator
double AWGN_generator();
// update the trellis at every new codeword decoding
void update_trellis (decoder_node_t ***Trellis_old, double *Y, int num_sections, int num_states, int num_outputs, int list_size, int num_inputs, int num_next_states);
// function to set communication pipes between processes
void set_pipe(int *** fd, int num_proc);//set number of pipes to use
// read the generator matrix
void read_G_matrix (char *G_filename, unsigned int ***G, int k, int n);
// functions to free dynamical arrays
void free2D(int **Matrix, int num_rows);
void free2D_unsigned(unsigned int **Matrix, int num_rows);
void free_Trellis (decoder_node_t **Trellis_nodes, int num_sections, int num_states, int num_next_states, int list_size);

// ###########################################################
// --- MAIN FUNCTION
int main(int argc, char *argv[])
{
    
    char *trellis_filename, *H_filename, *output_filename, *buf, *gen_matrix_filename; // buf is used to communicate between processes

    int k, n, m, num_errors_per_process, count_errors, num_WAVA_iterations, list_size, num_sections, num_states, num_outputs, num_inputs, num_next_states, num_errors, *ML_message, **H, r, i_check, i_codeword, i_bit, stop, count_punct, punct_every, num_punctures, num_processes, i_process, **f_pipe, **q_pipe, num_codewords_per_process, status,i_u_bit, total_bits, i_row; // k is the number of info bits, n the number of codeword bits of the convolutional code without CRC, m is the degree of the CRC code, ML_message is the found maximum likelihood message (if anything is found), r is the number of rows of the parity check matrix, i_ are indexes used, stop is a flag, f_pipe and q_pipe are used by the processes to communicate, status is used to check the status of the processes
	
	unsigned int *u, *v, **G, *all_zero; // u contains the info bit sent (message), v is the codeword sent, G the generator matrix of the code, all_zero is a all_zero vector

    decoder_node_t **Trellis_nodes; // Trellis structure

    double *Y, start_EbN0, end_EbN0, step_EbN0, CER, i_EbN0, std_noise, temp_double_value; // Y is the modulated codeword, all EbN0 terms are in dB, CER is the Codeword Error Rate, std_noise the standard deviation of the noise, while temp_double_value is used as temporary variable

    FILE *file_out; // file pointer
    pid_t pid; // process ID for the multi-process
    
    // CHECK THAT THERE ARE ALL ARGUMENTS AND SAVE THEM, OTHERWISE RETURN AN ERROR
    if (argc != 16){
        printf("ERROR! The program arguments shall be:\n./PLVA trellis_filename H_filename output_filename num_WAVA_iterations list_size num_errors k n m punctures num_processes G_matrix_filename start_EbN0 end_EbN0 step_EbN0\n");
        exit(1);
    }
    else {
        trellis_filename=argv[1];
        H_filename=argv[2];
        output_filename=argv[3];
        num_WAVA_iterations=atoi(argv[4]); // if 1, pure Viterbi
        list_size=atoi(argv[5]);
        num_errors=atoi(argv[6]); // number of errors to be counted at each E_b/N_0 simulated point
        k=atoi(argv[7]);
        n=atoi(argv[8]); 
        m=atoi(argv[9]); // degree of the CRC code
        num_punctures=atoi(argv[10]);
        num_sections=(k+m);
        num_outputs=(n/k); // number of output symbols
	    num_processes=atoi(argv[11]);
		gen_matrix_filename = argv[12];
        start_EbN0=atof(argv[13]);
        end_EbN0=atof(argv[14]);
        step_EbN0=atof(argv[15]);
    }
    
    // allocate memory and check if there are any errors
    Y = (double*)malloc(num_sections * num_outputs * sizeof(double));
	v = (unsigned int*)malloc(num_sections * num_outputs * sizeof(unsigned int));
	all_zero = (unsigned int*)malloc(num_sections * num_outputs * sizeof(unsigned int));
    u = (unsigned int*)malloc(k * sizeof(unsigned int));
    if(Y==NULL || u == NULL || v == NULL || all_zero == NULL){
        printf("ERROR! While allocating memory\n");
        exit(2);
    }
	
	// read the generator matrix from file
 	read_G_matrix (gen_matrix_filename, &G, k, num_sections * num_outputs);

    // set puncturing period
    punct_every = (int)floor((double)(num_sections) / (double)num_punctures);
    
    // initialize pipes for communications among processes
    set_pipe(&f_pipe, num_processes);//on f_pipe children wait
    set_pipe(&q_pipe, num_processes);//on q_pipe children send results to the parent

    // set number of errors to be computed by each process
    num_errors_per_process = (int)ceil(num_errors/num_processes); 

    // initialize the all-zero vector and the vector with the modulated symbols
	for (i_bit = 0; i_bit < num_sections * num_outputs; i_bit++) all_zero[i_bit] = (unsigned int) 0;
    for(i_bit=0;i_bit <num_sections*num_outputs;i_bit++) Y[i_bit]= (double)0.0;

    // build the trellis structure
    Trellis_nodes = initialize_trellis (trellis_filename, Y, num_sections, &num_states, num_outputs, 2, &num_inputs, &num_next_states) ;
    
    // allocate memory for the maximum likelihood message
    ML_message=(int *)malloc(num_inputs*num_sections*sizeof(int));
    if(ML_message==NULL){
        printf("ERROR! While allocating memory\n");
        exit(2);
    }

    // copy the parity check matrix
    read_H(H_filename, &H, num_inputs*num_sections, &r);
	
	total_bits = num_sections*num_outputs;
    
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
                
                // for loop over the different values of E_b / N_0
                for (i_EbN0 = start_EbN0; i_EbN0 <= end_EbN0; i_EbN0 += step_EbN0){
                    // initialize variables and compute noise standard deviation
                    CER = 0.0;
                    std_noise=sqrt(1.0/(2.0*((double)k/(double)((k+m-num_punctures)*num_outputs/num_inputs))*db2dec(i_EbN0)));
                    count_errors = 0;
                    i_codeword = 0;
					
					for (i_row = 0; i_row < k; i_row++) u[i_row]=0;
					
                    // while loop of encoding -> channel -> decoding until num_errors_per_process errors are counted by the process
            		while (count_errors < num_errors_per_process){
						
                        // generate random message u and encode it via G to obtain the transmitted codeword v
						memcpy(v,all_zero,total_bits*sizeof(unsigned int));
						for (i_row = 0; i_row < k; i_row++){
							u[i_row] = (unsigned int) rand() % 2;
							if (u[i_row] == 1){
								for (i_bit = 0; i_bit < total_bits; i_bit++) v[i_bit]+=G[i_row][i_bit];
							}							
						}
						for (i_bit = 0; i_bit < total_bits; i_bit++) v[i_bit]=v[i_bit]%2;
						
                        // generate transmitted modulated codeword via BSPK mapping (0 -> +1, 1 -> -1), puncture it (0->0,1->0) and add noise
                        count_punct=0;
                        for(i_bit=0;i_bit<total_bits;i_bit++){
                            if (i_bit!=0 && num_punctures>0 && ((i_bit%punct_every) == 0) && (count_punct < (num_punctures*num_outputs/num_inputs))){
                                Y[i_bit]= (double)0.0;
                                count_punct++;
                            }
                            else (v[i_bit] == 0) ? (Y[i_bit]= 1.0 + (std_noise * AWGN_generator())) : (Y[i_bit]= -1.0 + (std_noise * AWGN_generator()));
                        }
                        
                        // update the weigths of the edges of the trellis with the distance from the received codeword
                        update_trellis (&Trellis_nodes, Y, num_sections, num_states, num_outputs, 2, num_inputs, num_next_states) ;
                        
                        // if more than one WAVA iteration, run WAVA with list 1 to update the initial distances of the nodes in the initial trellis section. At the last WAVA iteration build the tree-trellis
                        WAVA_Update_Distances_and_Build_Tree_Trellis(&Trellis_nodes, num_sections, num_states, num_next_states, num_WAVA_iterations);
                            
                        // run SLVA
                        ML_message=Serial_List_Viterbi_Decoder(&Trellis_nodes, num_sections, num_states, num_next_states, list_size, num_inputs, H, r);
                        
                        // check if the found message is correct, otherwise increment the counter of the errors
                        stop=0;
                        for(i_bit=0;i_bit<num_sections*num_inputs-m && stop==0; i_bit++){
							if (ML_message[i_bit] != u[i_bit]){
								stop=1;
								count_errors++;
							}
                        }
                        
                        // free previously allocated memory
						free(ML_message);
						
                        // increment number of transmitted codewords
            			i_codeword++;

                    } // end while loop over the counter number of erroneus received codewords

                    // send to the parent process the number of transmitted codewords to reach num_errors_per_process errors
                    if(write(q_pipe[i_process][1],&i_codeword,sizeof(int))!=sizeof(int)){
                        fprintf(stderr,"ERROR while writing on pipe\n");
                        exit(99);
					}
                } // end for loops over E_b/N_0
                
                // close pipes to communicate with the parent and kill the child process
                close(q_pipe[i_process][1]);
                
                // free allocated memory
                free2D_unsigned(G, k);
                free2D(H, r);
                free2D(f_pipe,num_processes);
                free2D(q_pipe,num_processes);
                free(u);
                free(v);
                free(all_zero);
                free(Y);
                free_Trellis(Trellis_nodes, num_sections, num_states, num_next_states, 2);
                
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
        // collect the result of each process and sum them
        for(i_process=0; i_process<num_processes; i_process++){
            read(q_pipe[i_process][0],&num_codewords_per_process,sizeof(int));
            CER = CER + (double)num_codewords_per_process;
        }
        CER = (double)(num_processes*num_errors_per_process) / CER;
        
        // print result for each E_b / N_0 in the file and in the standard output
        fprintf(file_out,"%f\n",CER);
        printf("Eb/N0 = %.2f , CER = %f\n",i_EbN0, CER);
    }

    // when processed all E_b / N_0 close the output file and wait that each child process is terminated correctly
    fclose(file_out);
    for(i_process=0; i_process<num_processes; i_process++){//read on q_pipe
        while(read(q_pipe[i_process][0],&temp_double_value,sizeof(double))>0){}
        close(q_pipe[i_process][0]);
        wait(&status);
    }

    // free allocated memory
    //free(ML_message);
	free2D_unsigned(G, k);
    free2D(H, r);
    free2D(f_pipe,num_processes);
    free2D(q_pipe,num_processes);
	free(u);
	free(v);
	free(all_zero);
    free(Y);
    free_Trellis(Trellis_nodes, num_sections, num_states, num_next_states, 2);
    
    // return that all is ok
    return 0;
}

// ##################################
// --- functions description

// function to set the trellis structure reading from a file the structure
decoder_node_t **initialize_trellis (char *trellis_filename, double *Y, int num_sections, int *num_states, int num_outputs, int list_size, int *num_inputs, int *num_next_states){
    /*
    The file containing the Trellis must be of the form :
        num_states num_inputs num_next_states
        next_state0 next_state1 ... previous_state0 previous_state1 ... output_state0_0 output_state0_1 ...  input_state0_0 input_state0_1 ... // i_state = 0
        next_state0 next_state1 ... previous_state0 previous_state1 ... output_state0_0 output_state0_1 ...  input_state0_0 input_state0_1 ... // i_state = 1
        ...
        next_state0 next_state1 ... previous_state0 previous_state1 ... output_state0_0 output_state0_1 ...  input_state0_0 input_state0_1 ... // i_state = num_states -1
    */
    
    // define variables
    FILE *file_trellis;
    decoder_node_t **Trellis_nodes;
    int i_section, i_state, i, j, n_states, n_inputs, n_next_states;
    double temp_edge_cost, temp_value;
    
    // open file   
    file_trellis = fopen(trellis_filename, "r");
    if (file_trellis) {
        // read file info about trellis structure
        fscanf(file_trellis, "%d %d %d", &n_states, &n_inputs, &n_next_states);
        
        // allocate memory for the trellis and check that allocation is ok
        Trellis_nodes = (decoder_node_t **)malloc((num_sections+1)*sizeof(decoder_node_t*));
        if (Trellis_nodes == NULL) {
            printf("ERROR! Problems in allocating memory to the Trellis structure\n");
            exit(5);
        }
        
        // for each trellis section set connections, weights, inputs, ecc...
        for (i_section=0 ; i_section < num_sections+1 ; i_section++){
            Trellis_nodes[i_section] = (decoder_node_t *)malloc(n_states*sizeof(decoder_node_t));
            if (Trellis_nodes == NULL) {
                printf("ERROR! Problems in allocating memory to the Trellis structure\n");
                exit(6);
            }

            // for each trellis state of a section set connections, weights, inputs, ecc...
            for (i_state = 0; i_state < n_states ; i_state++){
                Trellis_nodes[i_section][i_state].next_states = (int*) malloc (n_next_states*sizeof(int));
                Trellis_nodes[i_section][i_state].outputs = (double**) malloc (n_next_states*sizeof(double*));
                Trellis_nodes[i_section][i_state].inputs = (int**) malloc (n_next_states*sizeof(int*));
                Trellis_nodes[i_section][i_state].edge_costs = (double*) malloc (n_next_states*sizeof(double));
                
                if (Trellis_nodes[i_section][i_state].next_states != NULL || Trellis_nodes[i_section][i_state].outputs != NULL || Trellis_nodes[i_section][i_state].inputs != NULL || Trellis_nodes[i_section][i_state].edge_costs != NULL)
                {
                    // copy next states from i_state
                    for (i=0 ; i< n_next_states; i++){
                        if (i_section == 0) fscanf(file_trellis, "%d", &Trellis_nodes[i_section][i_state].next_states[i]);
                        else Trellis_nodes[i_section][i_state].next_states[i] = Trellis_nodes[i_section-1][i_state].next_states[i];
                    }
                    // copy outputs and compute edge costs
                    for (i=0 ; i< n_next_states; i++){
                        Trellis_nodes[i_section][i_state].outputs[i] = (double*) malloc (num_outputs*sizeof(double));
                        if (Trellis_nodes[i_section][i_state].outputs[i] != NULL){
                            temp_edge_cost = 0;
                            for (j=0; j<num_outputs ; j++){
                                if (i_section == 0) fscanf(file_trellis, "%lf", &Trellis_nodes[i_section][i_state].outputs[i][j]);
                                else Trellis_nodes[i_section][i_state].outputs[i][j] = Trellis_nodes[i_section-1][i_state].outputs[i][j];
                                
                                temp_value=Y[i_section*num_outputs+j] - Trellis_nodes[i_section][i_state].outputs[i][j];
                                temp_edge_cost = temp_edge_cost + temp_value*temp_value;
                            }
                            Trellis_nodes[i_section][i_state].edge_costs[i] = temp_edge_cost;
                        }
                        else{
							printf("ERROR! While allocating memory for the Trellis\n");
							exit(7);
						}
                    }
                    // copy inputs
                    for (i=0 ; i< n_next_states; i++){
                        Trellis_nodes[i_section][i_state].inputs[i] = (int*) malloc (n_inputs*sizeof(int));
                        if (Trellis_nodes[i_section][i_state].inputs[i] != NULL){
                            for (j=0; j<n_inputs ; j++){
                                if (i_section == 0) fscanf(file_trellis, "%d", &Trellis_nodes[i_section][i_state].inputs[i][j]);
                                else Trellis_nodes[i_section][i_state].inputs[i][j]=Trellis_nodes[i_section-1][i_state].inputs[i][j];
                            }
                        }
                        else{
							printf("ERROR! While allocating memory for the Trellis\n");
							exit(8);
						}
                    }
                    // generate list infos
                    for (i=0 ; i<list_size; i++){
                        if (i_section == 0) {
                            Trellis_nodes[i_section][i_state].starting_node[i] = i_state;
                            Trellis_nodes[i_section][i_state].path_metric[i] = 0;
                        }
                        else {
                            Trellis_nodes[i_section][i_state].starting_node[i] = -1;
                            Trellis_nodes[i_section][i_state].path_metric[i] = -1;
                        }
                        Trellis_nodes[i_section][i_state].previous_node[i] = -1;
                        Trellis_nodes[i_section][i_state].previous_index[i] = -1;

                    }
                }
                else{
					printf("ERROR! While allocating memory for the Trellis\n");
					exit(9);
				}
                
            }
        }
            
        // close read file
        fclose(file_trellis);

        // set read parameters
        *num_inputs=n_inputs;
        *num_next_states=n_next_states;
        *num_states=n_states;
        
    }
    else{ // if problem in opening the file
        printf("ERROR! Problems in reading the file : <%s>\n", trellis_filename);
        exit(10);
    }
    
    return Trellis_nodes;
}

// update the weights of the edges of the trellis structure with the one of the new received vector
void update_trellis (decoder_node_t ***Trellis_old, double *Y, int num_sections, int num_states, int num_outputs, int list_size, int num_inputs, int num_next_states){

    // define variables
    decoder_node_t **Trellis_nodes = * Trellis_old;
    int i_section, i_state, i, j;
    double temp_edge_cost, temp_value;
    
    // compute the new weights of each edge and store them in the structure
    for (i_section = 0; i_section < num_sections; i_section ++){
        for(i_state = 0; i_state < num_states; i_state++){
            for (i=0 ; i< num_next_states; i++){
                temp_edge_cost = 0;
                for (j=0; j< num_outputs; j++){
                    temp_value=Y[i_section*num_outputs+j] - Trellis_nodes[i_section][i_state].outputs[i][j];
                    temp_edge_cost += temp_value*temp_value;
                }
                Trellis_nodes[i_section][i_state].edge_costs[i] = temp_edge_cost;
            }
            // update list infos
            if (i_section == 0) {
                Trellis_nodes[0][i_state].starting_node[0] = i_state;
                Trellis_nodes[0][i_state].path_metric[0] = 0.0;
                Trellis_nodes[0][i_state].starting_node[1] = i_state;
                Trellis_nodes[0][i_state].path_metric[1] = 0.0;
            }
            else {
                Trellis_nodes[i_section][i_state].starting_node[0] = -1;
                Trellis_nodes[i_section][i_state].path_metric[0] = -1;
                Trellis_nodes[i_section][i_state].starting_node[1] = -1;
                Trellis_nodes[i_section][i_state].path_metric[1] = -1;
            }
            Trellis_nodes[i_section][i_state].previous_node[0] = -1;
            Trellis_nodes[i_section][i_state].previous_index[0] = -1;
            Trellis_nodes[i_section][i_state].previous_node[1] = -1;
            Trellis_nodes[i_section][i_state].previous_index[1] = -1;
        }
    }
    *Trellis_old = Trellis_nodes;
}


/*
wrap-around Viterbi decoder for the first (num_WAVA - 1) iterations of the WAVA algorithm, to set the initial distances of the nodes in first section
plus this function builds the tree-trellis structure: each node has a list of 2 made of the best incoming path in the first entry, plus path from the other incoming edge in the second entry
*/
void WAVA_Update_Distances_and_Build_Tree_Trellis(decoder_node_t ***New_Trellis_nodes, int num_sections, int num_states, int num_next_states, int num_WAVA){
    // declare variables
    int i_section, i_state, j_state, index, i_WAVA, i_value[num_states], min_index, other_index;
    double values[num_states][num_next_states], other_path_metric;
    decoder_node_t **Trellis_nodes = *New_Trellis_nodes;
	
    // set to 0 the path metric at each node in the last trellis section
    for (i_state = 0; i_state < num_states ; i_state++) Trellis_nodes[num_sections][i_state].path_metric[0]=0.0;
    
    // for loop over the WAVA iterations
    for (i_WAVA=0 ; i_WAVA < num_WAVA ; i_WAVA++){
        // update every state in the first trellis section at every new WAVA iteration, with the values of the same nodes in the final trellis section 
        for (i_state = 0; i_state < num_states ; i_state++){
            Trellis_nodes[0][i_state].starting_node[0] = i_state;
            Trellis_nodes[0][i_state].previous_node[0] = -1;
            Trellis_nodes[0][i_state].previous_index[0] = -1;
            Trellis_nodes[0][i_state].path_metric[0] = Trellis_nodes[num_sections][i_state].path_metric[0];
        }
        
		// run Viterbi and store both best path and other incoming path to each node at each trellis section 
        for (i_section=0 ; i_section < num_sections ; i_section++){			
				for (i_state = 0; i_state < num_states ; i_state++) i_value[i_state] = 0;
				
				for (i_state = 0; i_state < num_states ; i_state++){
					for (j_state = 0; j_state < num_next_states; j_state++){
                        // index of a node connected to i_state node
						index = Trellis_nodes[i_section][i_state].next_states[j_state];
						
                        // check if it is the first time you enter in this node
						if (i_value[index] == 0){
							Trellis_nodes[i_section+1][index].starting_node[0] = Trellis_nodes[i_section][i_state].starting_node[0];
							Trellis_nodes[i_section+1][index].previous_node[0] = i_state;
							Trellis_nodes[i_section+1][index].previous_index[0] = j_state;
							Trellis_nodes[i_section+1][index].path_metric[0] = Trellis_nodes[i_section][i_state].path_metric[0] + Trellis_nodes[i_section][i_state].edge_costs[j_state];
						}
						else {
							other_path_metric = Trellis_nodes[i_section][i_state].path_metric[0] + Trellis_nodes[i_section][i_state].edge_costs[j_state];
							
							// if second incoming distance is bigger, store this path data in the second position of the list 
							if (Trellis_nodes[i_section+1][index].path_metric[0] < other_path_metric){
								
								Trellis_nodes[i_section+1][index].starting_node[1] = Trellis_nodes[i_section][i_state].starting_node[0];
								Trellis_nodes[i_section+1][index].previous_node[1] = i_state;
								Trellis_nodes[i_section+1][index].previous_index[1] = j_state;
								Trellis_nodes[i_section+1][index].path_metric[1] = other_path_metric;
							}
							// if second incoming distance is shorter than the one saved before, store this path data in the first position of the list and move the previous in the second position
							else{
								// move previous stored
								Trellis_nodes[i_section+1][index].starting_node[1] = Trellis_nodes[i_section+1][index].starting_node[0];
								Trellis_nodes[i_section+1][index].previous_node[1] = Trellis_nodes[i_section+1][index].previous_node[0];
								Trellis_nodes[i_section+1][index].previous_index[1] = Trellis_nodes[i_section+1][index].previous_index[0];
								Trellis_nodes[i_section+1][index].path_metric[1] = Trellis_nodes[i_section+1][index].path_metric[0];
								// copy new								
								Trellis_nodes[i_section+1][index].starting_node[0] = Trellis_nodes[i_section][i_state].starting_node[0];
								Trellis_nodes[i_section+1][index].previous_node[0] = i_state;
								Trellis_nodes[i_section+1][index].previous_index[0] = j_state;
								Trellis_nodes[i_section+1][index].path_metric[0] = other_path_metric;
							}
								
						}
						
						i_value[index]++;
					}// close neighboors
				} // close state for loop
        } // close section for loop
        
		
    } // close num_WAVA for loop
    
    *New_Trellis_nodes = Trellis_nodes; 

}

/*
Serial List Viterbi Algorithm run at the last WAVA iteration
*/
int *Serial_List_Viterbi_Decoder(decoder_node_t ***Trellis_New_nodes, int num_sections, int num_states, int num_next_states, int list_size, int num_inputs, int **H, int r){
    // define variables
    int i_section, i_state, i_list, i_list_state, actual_list_size, next_list_size, stop = 0, best_codeword_index, best_state, *temp_codeword, *best_codeword, counter, j_state, index, last_section, found_a_codeword, i_stack;
    decoder_node_t **Trellis_nodes = *Trellis_New_nodes;
	
	heap_t *heap_stack=(heap_t *)calloc(1, sizeof (heap_t));

	stack_node_t stack_node, best_node;
    int stack_size=0,temp_stack1_size=0, temp_stack2_size=0, list_index=0, new_split, starting_section, starting_state, count_list ;
	double worst_path_metric, temp_path_metric, partial_path_metric;

    found_a_codeword = 0;
    actual_list_size = 1;
    next_list_size = 1;
		
	// store in the heap datastructure all the paths found at the nodes in the last trellis section
    // remember that each node in the stack represents a path and it has the following data:
    // (starting node, ending node, path metric, number of splits from the best path, section indexes of the splits, nodes of the splits, partial path metric from the last split to the ending node)
	stack_size=0;
	worst_path_metric = -1.0;
    for (i_state = 0; i_state < num_states ; i_state++){
		stack_node.starting_node = Trellis_nodes[num_sections][i_state].starting_node[0];
		stack_node.ending_node = i_state;
		stack_node.num_splits = 0;
		stack_node.path_metric = Trellis_nodes[num_sections][i_state].path_metric[0];
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
		temp_codeword = back_Trellis_stack(Trellis_nodes,best_node,num_sections,num_next_states, num_inputs);
		//check CRC condition via syndrome check
		if (check_syndrome(temp_codeword,H,num_inputs*num_sections,r)==1){
            // save as best message
			found_a_codeword = 1;
			best_codeword = temp_codeword;
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
	while (found_a_codeword == 0 && count_list <= list_size && stack_size>0){
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
			starting_section = num_sections;
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
			temp_codeword = back_Trellis_stack(Trellis_nodes,best_node,num_sections,num_next_states, num_inputs);
			// check CRC syndrome
			if (check_syndrome(temp_codeword,H,num_inputs*num_sections,r)==1){
				found_a_codeword = 1;
				best_codeword= temp_codeword;
				
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
    if (found_a_codeword == 0){
	    best_codeword=(int*)malloc(num_sections*num_inputs*sizeof(int));
	    if(best_codeword==NULL){
            printf("ERROR! While allocating memory to the final codeword\n");
            exit(11);
        }
	    for (i_section = 0; i_section < num_sections*num_inputs ; i_section++) best_codeword[i_section] = -1;
    }

    *Trellis_New_nodes = Trellis_nodes;

    return best_codeword;
    
}

// go bacward in the trellis and extract the input message associated to a path stored in the stack
int *back_Trellis_stack(decoder_node_t **Trellis_nodes,stack_node_t stack_node,int num_sections,int num_next_states, int num_inputs){

    int i_section, i_input, i_next_state, counter_splits, initial_state, previous_state, actual_state, previous_list_index, i_bit, *codeword, stop=0;

    codeword=(int *)malloc(num_sections*num_inputs*sizeof(int));
    if(codeword == NULL){
        printf("ERROR while allocating space for the codeword!");
        exit(12);
    }

    initial_state = stack_node.ending_node;
    actual_state = initial_state;

    i_bit=num_sections - 1;
	
	
	if (stack_node.num_splits == 0) {
		for (i_section = num_sections; i_section > 0 ; i_section--){
			codeword[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[0];
			i_bit --;
			actual_state = Trellis_nodes[i_section][actual_state].previous_node[0];
		}
	}
	else {
		counter_splits = 0;
		for (i_section = num_sections; i_section > 0 ; i_section--){
			if (i_section == stack_node.section_split_indexes[counter_splits]){
				codeword[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[1];
				i_bit --;
				actual_state = Trellis_nodes[i_section][actual_state].previous_node[1];
				if (counter_splits  < stack_node.num_splits) counter_splits++;
			}
			else{
				codeword[i_bit]=Trellis_nodes[i_section][actual_state].previous_index[0];
				i_bit --;
				actual_state = Trellis_nodes[i_section][actual_state].previous_node[0];
			}
		}
		
	}
	
	
	
	return codeword;
}

int check_syndrome(int *temp_codeword, int **H, int n, int r){
    int syndrome[r], i_bit, i_row, stop=0;
    
    for(i_row=0; i_row<r && stop==0; i_row++){
        syndrome[i_row]=0;
        for(i_bit = 0; i_bit < n; i_bit++) syndrome[i_row] = syndrome[i_row]+H[i_row][i_bit]*temp_codeword[i_bit];
        syndrome[i_row] = syndrome[i_row]%2;
        if (syndrome[i_row] != 0) stop=1;
    }
    
    return 1-stop;
    /* return 1 if syndrome-check is passed with allzero syndrome vector,
    return 0 if syndrome vector is not the allzero vector*/
}

// function to read the parity check matrix of the outer code from a file
void read_H(char *H_filename, int ***H_main, int n, int *r){
/*
    The file containing Y must be of the form :
    
        num_rows
        H_00 H_01 ... H_0n
        ...
        ... H_rn
    
    where n = num_columns, r = num_rows
    */
    FILE *file_H;
    int i_row, i_column, num_rows, **H;
    
    file_H = fopen(H_filename, "r");
    if (file_H) {
        fscanf(file_H, "%d", &num_rows);
        
        H = (int **)malloc(num_rows*sizeof(int*));
        if (H == NULL) {
            printf("ERROR! Problems in allocating memory to H");
            exit(13);
        }
        else{
            for(i_row=0;i_row<num_rows; i_row++){
                H[i_row]=(int*)malloc(n*sizeof(int));
                if (H[i_row] == NULL) {
                    printf("ERROR! Problems in allocating memory to H");
                    exit(14);
                }
                else{
                    for (i_column=0; i_column<n; i_column++){
                        fscanf(file_H, "%d",&H[i_row][i_column]);
                    }
                }
            }
        }        
        
        fclose(file_H);
    }
    else{
        printf("ERROR! Problems in reading the file : <%s>\n", H_filename);
        exit(15);
    }
    
    *H_main = H;
    *r = num_rows;
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

// function to read from file the generator matrix
void read_G_matrix (char *G_filename, unsigned int ***G_pointer, int k, int n){
/*
    The file containing G must be of the form :
    
        G_0,0 G_0,1 ... G_0,n
        ...
        ... G_k-1,n
    
    */
    FILE *file_G;
    int i_row, i_column;
	unsigned int **G;
    
    file_G = fopen(G_filename, "r");
    if (file_G) {
        
        G = (unsigned int **)malloc(k*sizeof(unsigned int*));
        if (G == NULL) {
            printf("ERROR! Problems in allocating memory to G");
            exit(18373);
        }
        else{
            for(i_row=0;i_row<k; i_row++){
                G[i_row]=(unsigned int*)malloc(n*sizeof(unsigned int));
                if (G[i_row] == NULL) {
                    printf("ERROR! Problems in allocating memory to G");
                    exit(18464);
                }
                else{
                    for (i_column=0; i_column<n; i_column++){
                        fscanf(file_G, "%d",&G[i_row][i_column]);
                    }
                }
            }
        }        
        
        fclose(file_G);
    }
    else{
        printf("ERROR! Problems in reading the file : <%s>\n", G_filename);
        exit(15364);
    }
    
    *G_pointer = G;
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
    int i, j, k;
    
	
    stack_node_t best_node = h->nodes[1];
	
    h->nodes[1] = h->nodes[h->len];
 
    h->len--;
 
    i = 1;
    while (i!=h->len+1) {
        k = h->len+1;
        j = 2 * i;
        if (j <= h->len && h->nodes[j].path_metric < h->nodes[k].path_metric) {
            k = j;
        }
        if (j + 1 <= h->len && h->nodes[j + 1].path_metric < h->nodes[k].path_metric) {
            k = j + 1;
        }
        h->nodes[i] = h->nodes[k];
        i = k;
    }
	
    return best_node;
}

// functions to deallocate memory
void free2D(int **Matrix, int num_rows){
    int i_row;
    
    for(i_row=0; i_row<num_rows; i_row++) free(Matrix[i_row]);
    free(Matrix);
}
void free_Trellis (decoder_node_t **Trellis_nodes, int num_sections, int num_states, int num_next_states, int list_size){
    int i_section,i_state,i,j;

    for(i_section=0; i_section<num_sections; i_section++){
        for(i_state=0; i_state<num_states; i_state++){
            free(Trellis_nodes[i_section][i_state].next_states);
            for(i=0; i<num_next_states; i++){
                free(Trellis_nodes[i_section][i_state].outputs[i]);
                free(Trellis_nodes[i_section][i_state].inputs[i]);
            }
            free(Trellis_nodes[i_section][i_state].edge_costs);
        }
        free(Trellis_nodes[i_section]);
    }
    free(Trellis_nodes);
}
void free2D_unsigned(unsigned int **Matrix, int num_rows){
    int i_row;
    
    for(i_row=0; i_row<num_rows; i_row++) free(Matrix[i_row]);
    free(Matrix);
}
