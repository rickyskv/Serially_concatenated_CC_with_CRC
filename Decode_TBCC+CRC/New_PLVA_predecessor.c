// Parrallel List (Wrap-Around) Viterbi Algorithm

// Note: this file contains a simple implementation of the P-LVA.
// It is not multi-threading. The multiple processes are used to speed-up simulations, where each process is an indipendent decoder/simulator.
// A central processor receives from the other processes the number of codewords simulated before reaching a number of errors defined by the user.
// To set properly the parameters of this C code, use the associated MATLAB file.

// Possible improvements will be published in the future:
// - implement the sorter network proposed in the thesis to see if the compiler use parallel structures of the processor
// - quantize the decoder
// - use intel vector instructions to speed-up the simulator when fixed-point arithmetic is used
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
//#include <stdint.h>
#include <stdlib.h>
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
// Viterbi Algorithm
void Viterbi_Decoder(connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK);
// Parallel List Viterbi Algorithm
void Parallel_List_Viterbi_Decoder(connections_node_t Trellis_links[NumStates], decoder_node_t ***Trellis_nodes, int16_t list_size, uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK);
// backward function to find the codeword/message corresponding to a certain path
void back_Trellis_small(connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates],int i_state, int i_list, uint8_t (*codeword)[NumSections*NumInputs]);
// backward function to find the codeword/message corresponding to a certain list path
void back_Trellis(connections_node_t Trellis_links[NumStates], decoder_node_t ***Trellis_nodes,int i_state, int i_list, uint8_t (*codeword)[NumSections*NumInputs]);//void back_Trellis(connections_node_t Trellis_links[NumStates], decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], decoder_list_t list_entry, uint8_t (*codeword)[NumSections*NumInputs]);
// Additive-White Gaussian Noise generator
double AWGN_generator();
// update the trellis at every new codeword decoding
void update_trellis (connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], double Y[NumSections*NumOutputs], int16_t list_size);
// function to set communication pipes between processes
void set_pipe(int *** fd, int num_proc);//set number of pipes to use
// read a puncturing pattern
void read_puncturing(char *Puncturing_file, uint8_t (*puncturing_v)[NumSections*NumOutputs]);
// functions to free dynamical arrays
void free2D(int **Matrix, int num_rows);
void free2D_unsigned(unsigned int **Matrix, int num_rows);
void free2D_Trellis(decoder_node_t **Trellis_nodes);

// ###########################################################
// --- MAIN FUNCTION
int main(int argc, char *argv[])
{
    
    char *output_filename, *buf, *puncturing_vector_filename; // buf is used to communicate between processes

    int num_errors_per_process, count_errors, count_UE, count_NACK, tot_NACK, num_errors, temp_list, i_list, i_section, i_state, i_WAVA, i_check, i_codeword, i_bit, stop, count_punct, punct_every, num_punctures, num_processes, i_process, **f_pipe, **q_pipe, num_codewords_per_process, status,i_u_bit, total_bits, i_row, total_list_size, actual_list_size;; // k is the number of info bits, n the number of codeword bits of the convolutional code without CRC, m is the degree of the CRC code, ML_message is the found maximum likelihood message (if anything is found), r is the number of rows of the parity check matrix, i_ are indexes used, stop is a flag, f_pipe and q_pipe are used by the processes to communicate, status is used to check the status of the processes
    uint8_t NACK, ML_message[NumSections * NumInputs], num_WAVA_iterations;
    int16_t list_size;
	
    uint8_t u[k], u_CRC[k+m], v[NumSections * NumOutputs], all_zero[NumSections * NumOutputs], puncturing_v[NumSections * NumOutputs]; // u contains the info bit sent (message), v is the codeword sent, all_zero is a all_zero vector
    
    decoder_node_t **Trellis_nodes;//[NumSections+1][NumStates]; // Trellis structure for PLVA (dinamically allocated in the heap for memory reasons)
    small_decoder_node_t Trellis_nodes_small[NumSections+1][NumStates]; // Trellis structure for VA/WAVA no list
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
    


    initialize_trellis(&Trellis_links);//, &Trellis_nodes);
    Trellis_nodes = (decoder_node_t**)malloc((NumSections+1)*sizeof(decoder_node_t*));
    if (Trellis_nodes==NULL){
        printf("ERROR while allocating memory!\n");
        exit(2);
    }
    for (i_section = 0; i_section <= NumSections; i_section++){
        Trellis_nodes[i_section] = (decoder_node_t*)malloc((NumStates)*sizeof(decoder_node_t));
        if (Trellis_nodes==NULL){
            printf("ERROR while allocating memory!\n");
            exit(3);
        }
    }

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
    for (i_bit = 0; i_bit < NumSections * NumOutputs; i_bit++) all_zero[i_bit] = (uint8_t) 0;
    for(i_bit=0;i_bit <NumSections*NumOutputs;i_bit++) Y[i_bit]= (double)0.0;
	
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
                
                for(i_state = 0; i_state < NumStates; i_state++){
                    Trellis_nodes[0][i_state].starting_node[0] = i_state;
                    for(i_list = 1; i_list < list_size; i_list++){
                        Trellis_nodes[0][i_state].starting_node[i_list] = -1;
                        Trellis_nodes[0][i_state].path_metric[i_list] = 10000.0;
                    }
                }
                temp_list=1;
                for(i_section = 1; i_section < (int)log2(list_size); i_section++){
                    temp_list<<=1;
                    for(i_state = 0; i_state < NumStates; i_state++){
                        for(i_list = temp_list; i_list < list_size; i_list++){
                            Trellis_nodes[i_section][i_state].starting_node[i_list] = -1;
                            Trellis_nodes[i_section][i_state].path_metric[i_list] = 10000.0;
                        }
                    }
                }

                // for loop over the different values of E_b / N_0
                for (i_EbN0 = start_EbN0; i_EbN0 <= end_EbN0; i_EbN0 += step_EbN0){
                    // initialize variables and compute noise standard deviation
                    CER = 0.0;
                    total_list_size = 0;
                    std_noise=sqrt(1.0/(2.0*((double)k/(((double)(k+m)*(double)NumOutputs/(double)NumInputs)-(double)num_punctures)*db2dec(i_EbN0))));
                    
                    count_errors = 0;
                    count_UE = 0;
                    count_NACK = 0;
                    i_codeword = 0;
					
                    // while loop of encoding -> channel -> decoding until num_errors_per_process errors are counted by the process
            	    while (count_errors < num_errors_per_process){
						
                        // generate random message u and encode it via G to obtain the transmitted codeword v
			for (i_bit=0;i_bit<k;i_bit++) u[i_bit]=((uint8_t)rand()) % 2;	
			CRC_encode (u, &u_CRC);
			TBCC_encode (u_CRC, &v);
			
                        // generate transmitted modulated codeword via BSPK mapping (0 -> +1, 1 -> -1), puncture it (0->0,1->0) and add noise
                        for(i_bit=0;i_bit<total_bits;i_bit++){
                            if (puncturing_v[i_bit]==0){
                            	Y[i_bit]= (double)0.0;
                            }
                            else {
                            	(v[i_bit] == 0) ? (Y[i_bit]= 1.0 + (std_noise * AWGN_generator())) : (Y[i_bit]= -1.0 + (std_noise * AWGN_generator()));
                            }
                        }
                        
                        // update the weigths of the edges of the trellis with the distance from the received codeword
                        update_trellis(Trellis_links, &Trellis_nodes_small, Y, 1) ;
                        
                        // if more than one WAVA iteration, run WAVA with list 1 to update the initial distances of the nodes in the initial trellis section

                        // the last WAVA iteration run PLVA
                        NACK=1;
                        for (i_WAVA = 0; i_WAVA<num_WAVA_iterations && NACK!=0;i_WAVA++){
                            if(i_WAVA<num_WAVA_iterations-1){
                                total_list_size+= 1;
                                Viterbi_Decoder(Trellis_links, &Trellis_nodes_small, &ML_message, &NACK);
                            }
                            else{
                                total_list_size+= list_size;
                                for(i_state=0;i_state<NumStates;i_state++){
                                    Trellis_nodes[NumSections][i_state].path_metric[0] = Trellis_nodes_small[NumSections][i_state].path_metric[0];
                                }
                                for(i_section=0;i_section<=NumSections;i_section++){
                                    for(i_state=0;i_state<NumStates;i_state++){
                                        Trellis_nodes[i_section][i_state].edge_costs[0] = Trellis_nodes_small[i_section][i_state].edge_costs[0];
                                        Trellis_nodes[i_section][i_state].edge_costs[1] = Trellis_nodes_small[i_section][i_state].edge_costs[1];
                                    }
                                }
	                        Parallel_List_Viterbi_Decoder(Trellis_links, &Trellis_nodes, list_size, &ML_message, &NACK);
	                    }
                        }
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
                        exit(99);}
                } // end for loops over E_b/N_0
                
                // close pipes to communicate with the parent and kill the child process
                close(q_pipe[i_process][1]);
                
                // free allocated memory
                
                free2D(f_pipe,num_processes);
                free2D(q_pipe,num_processes);
                free2D_Trellis(Trellis_nodes);
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
        total_list_size =0;
        // collect the result of each process and sum them
        for(i_process=0; i_process<num_processes; i_process++){
            read(q_pipe[i_process][0],&num_codewords_per_process,sizeof(int));
            CER = CER + (double)num_codewords_per_process;
            read(q_pipe[i_process][0],&actual_list_size,sizeof(int));
            total_list_size += actual_list_size;
            read(q_pipe[i_process][0],&count_NACK,sizeof(int));
            tot_NACK +=  count_NACK;
        }
        average_list_size = ((double)total_list_size)/CER;
        CER = (double)(num_processes*num_errors_per_process) / CER;
        
        // print result for each E_b / N_0 in the file and in the standard output
        //fprintf(file_out,"%f\n",CER);
        fprintf(file_out,"Eb/N0 = %.2f , CER = %.9f, E[L]= %.2f, NACK/(NACK+UE) %d/%d\n",i_EbN0, CER, average_list_size, tot_NACK, num_processes*num_errors_per_process);
        printf("Eb/N0 = %.2f , CER = %.9f, E[L]= %.2f, NACK/(NACK+UE) %d/%d\n",i_EbN0, CER, average_list_size, tot_NACK, num_processes*num_errors_per_process);
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
    free2D_Trellis(Trellis_nodes);
    
    // return that all is ok
    return 0;
}

// ##################################
// --- functions description

// update the weights of the edges of the trellis structure with the one of the new received vector
// update the weights of the edges of the trellis structure with the one of the new received vector
void update_trellis (connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], double Y[NumSections*NumOutputs], int16_t list_size){

    // define variables

    int16_t i_section, i_state, i_list, i, j;
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
    
    for (i_state = 0; i_state < NumStates ; i_state++){
        (*Trellis_nodes)[0][i_state].starting_node[0] = i_state;
        (*Trellis_nodes)[0][i_state].previous_node[0] = -1;
        (*Trellis_nodes)[0][i_state].previous_index[0] = -1;
        (*Trellis_nodes)[NumSections][i_state].path_metric[0]=0.0;
        (*Trellis_nodes)[0][i_state].path_metric[0]=0.0;
    }

}




/*
Parallel List Viterbi Algorithm run at the last WAVA iteration
*/
void Parallel_List_Viterbi_Decoder(connections_node_t Trellis_links[NumStates], decoder_node_t ***Trellis_nodes, int16_t list_size, uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK){

    // define variables
    int16_t i_bit, i_section, i_state, i_list, i_list_state, actual_list_size, next_list_size, temp_list, stop = 0, best_codeword_index, best_state, counter, j_state, last_section, found_a_codeword,found_TBCC = 0;
    uint8_t temp_codeword[NumSections*NumInputs], bit[NumStates][2];
	
    int16_t list_counter[NumStates][2], min_index[NumStates], index[NumStates][2];
	
    double best_path_metric, metric[NumStates][2];
    
    // initialize values
    found_a_codeword = 0;
    best_path_metric = -1;
    actual_list_size = 1;
    next_list_size = 1;

    for(i_state = 0; i_state < NumStates; i_state++){
	(*Trellis_nodes)[0][i_state].path_metric[0] = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
    }
    
    // list decoder
    // for loop over trellis sections
    for (i_section=1 ; i_section <= NumSections ; i_section++){

        // update list size L (1st section L=1, 2nd section L=2, 3rd section L=4, ... until reached the maximum list size defined by the user)
        actual_list_size = next_list_size;
        next_list_size = ((NumNextStates*actual_list_size) > list_size) ? list_size : NumNextStates*actual_list_size ;
        
        // look at the two predecessors, update their path metrics and choose the total best L paths

        for (i_state = 0; i_state < NumStates ; i_state++){
        
            index[i_state][0] = Trellis_links[i_state].previous_states[0];
	    index[i_state][1] = Trellis_links[i_state].previous_states[1];
	    bit[i_state][0] = Trellis_links[i_state].previous_states_inputs[0][0];
	    bit[i_state][1] = Trellis_links[i_state].previous_states_inputs[1][0];
	    
	    list_counter[i_state][0] = 0;
            list_counter[i_state][1] = 0;
            
            metric[i_state][0] = (*Trellis_nodes)[i_section-1][index[i_state][0]].path_metric[list_counter[i_state][0]] + (*Trellis_nodes)[i_section-1][index[i_state][0]].edge_costs[bit[i_state][0]];
	    metric[i_state][1] = (*Trellis_nodes)[i_section-1][index[i_state][1]].path_metric[list_counter[i_state][1]] + (*Trellis_nodes)[i_section-1][index[i_state][1]].edge_costs[bit[i_state][1]];

	    for (i_list = 0; i_list < next_list_size; i_list++){
	        min_index[i_state] = (metric[i_state][0] < metric[i_state][1]) ? 0 : 1;
	        
	        (*Trellis_nodes)[i_section][i_state].starting_node[i_list] = (*Trellis_nodes)[i_section-1][index[i_state][min_index[i_state]]].starting_node[list_counter[i_state][min_index[i_state]]];
		(*Trellis_nodes)[i_section][i_state].previous_node[i_list] = min_index[i_state];
		(*Trellis_nodes)[i_section][i_state].previous_index[i_list] = list_counter[i_state][min_index[i_state]];
		(*Trellis_nodes)[i_section][i_state].path_metric[i_list] = metric[i_state][min_index[i_state]];
		
		list_counter[i_state][min_index[i_state]]++;
		metric[i_state][min_index[i_state]] = (*Trellis_nodes)[i_section-1][index[i_state][min_index[i_state]]].path_metric[list_counter[i_state][min_index[i_state]]] + (*Trellis_nodes)[i_section-1][index[i_state][min_index[i_state]]].edge_costs[bit[i_state][min_index[i_state]]];
		
	    }
        }
            
            if (i_section == NumSections){
                for (i_state = 0; i_state < NumStates; i_state++){
                    i_list=0;
		    stop=0;
		    counter = 0;
                    for (i_list = 0; i_list < next_list_size && stop==0; i_list++){
                        (*Trellis_nodes)[NumSections][i_state].path_metric[i_list] -= (*Trellis_nodes)[0][i_state].path_metric[0];
                        /*if (list_size == 1 && ((*Trellis_nodes)[NumSections][i_state].path_metric[0]<best_path_metric || best_path_metric == -1)){
                            best_path_metric = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
                        }*/
		        // check the path(s) of a node if you have not found a tail-biting codeword respecting the CRC condition or if the actual path in the list has smaller distance than the found codeword which respects the CRC condition
                        if((*Trellis_nodes)[NumSections][i_state].path_metric[i_list] < best_path_metric || found_a_codeword == 0){
                            // check the tail-biting condition
                            if ((*Trellis_nodes)[NumSections][i_state].starting_node[i_list] == i_state) {
								
				// if tail-biting go backward and extract the corresponding message
                                back_Trellis(Trellis_links,Trellis_nodes,i_state,i_list,&temp_codeword);
    
                                // check that the message satisfies the CRC condition
                                if (CRC_check(temp_codeword)==1){
                                    for (i_bit = 0; i_bit < NumSections*NumInputs ; i_bit++) (*best_codeword)[i_bit] = temp_codeword[i_bit];
				    found_a_codeword=1;
                                    best_path_metric = (*Trellis_nodes)[NumSections][i_state].path_metric[i_list];
    
                                    stop=1;
                                }
                            }
                        }
                        else stop=1;
                    }
                }
                        
            }
    }
    
    (found_a_codeword == 0) ? ((*NACK) = 1) : ((*NACK) = 0);
    
}

/*
Viterbi Algorithm
*/
void Viterbi_Decoder(connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates], uint8_t (*best_codeword)[NumSections*NumInputs], uint8_t *NACK){

    // define variables
    int16_t i_bit, i_section, i_state, i_list, i_list_state, actual_list_size, next_list_size, temp_list, stop = 0, best_codeword_index, best_state, counter, j_state, last_section, found_a_codeword,found_TBCC = 0;
    uint8_t temp_codeword[NumSections*NumInputs], bit[NumStates][2];
	
    int16_t min_index[NumStates], index[NumStates][2];
	
    double best_path_metric, metric[NumStates][2];
    
    // initialize values
    found_a_codeword = 0;
    best_path_metric = -1;
    actual_list_size = 1;
    next_list_size = 1;

    for(i_state = 0; i_state < NumStates; i_state++){
	(*Trellis_nodes)[0][i_state].path_metric[0] = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
    }
    
    // list decoder
    // for loop over trellis sections
    for (i_section=1 ; i_section <= NumSections ; i_section++){

        // look at the two predecessors, update their path metrics and choose the total best L paths

        for (i_state = 0; i_state < NumStates ; i_state++){
        
            index[i_state][0] = Trellis_links[i_state].previous_states[0];
	    index[i_state][1] = Trellis_links[i_state].previous_states[1];
	    bit[i_state][0] = Trellis_links[i_state].previous_states_inputs[0][0];
	    bit[i_state][1] = Trellis_links[i_state].previous_states_inputs[1][0];
            
            metric[i_state][0] = (*Trellis_nodes)[i_section-1][index[i_state][0]].path_metric[0] + (*Trellis_nodes)[i_section-1][index[i_state][0]].edge_costs[bit[i_state][0]];
	    metric[i_state][1] = (*Trellis_nodes)[i_section-1][index[i_state][1]].path_metric[0] + (*Trellis_nodes)[i_section-1][index[i_state][1]].edge_costs[bit[i_state][1]];


	    min_index[i_state] = (metric[i_state][0] < metric[i_state][1]) ? 0 : 1;
	        
	    (*Trellis_nodes)[i_section][i_state].starting_node[0] = (*Trellis_nodes)[i_section-1][index[i_state][min_index[i_state]]].starting_node[0];
	    (*Trellis_nodes)[i_section][i_state].previous_node[0] = min_index[i_state];
	    (*Trellis_nodes)[i_section][i_state].previous_index[0] = 0;
	    (*Trellis_nodes)[i_section][i_state].path_metric[0] = metric[i_state][min_index[i_state]];

        }
            
        if (i_section == NumSections){
            for (i_state = 0; i_state < NumStates; i_state++){
                (*Trellis_nodes)[NumSections][i_state].path_metric[i_list] -= (*Trellis_nodes)[0][i_state].path_metric[0];
                /*if (((*Trellis_nodes)[NumSections][i_state].path_metric[0]<best_path_metric || best_path_metric == -1)){
                    best_path_metric = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
                }*/
	        // check the path(s) of a node if you have not found a tail-biting codeword respecting the CRC condition or if the actual path in the list has smaller distance than the found codeword which respects the CRC condition
                if((*Trellis_nodes)[NumSections][i_state].path_metric[0] < best_path_metric || found_a_codeword == 0){
                    // check the tail-biting condition
                    if ((*Trellis_nodes)[NumSections][i_state].starting_node[0] == i_state) {
							
			// if tail-biting go backward and extract the corresponding message
                        back_Trellis_small(Trellis_links,Trellis_nodes,i_state,0,&temp_codeword);

                        // check that the message satisfies the CRC condition
                        if (CRC_check(temp_codeword)==1){
                            for (i_bit = 0; i_bit < NumSections*NumInputs ; i_bit++) (*best_codeword)[i_bit] = temp_codeword[i_bit];
			    found_a_codeword=1;
                            best_path_metric = (*Trellis_nodes)[NumSections][i_state].path_metric[0];
                        }
                    }
                }
            }
                        
        }
    }
    
    (found_a_codeword == 0) ? ((*NACK) = 1) : ((*NACK) = 0);
    
}

// go backward in the trellis and extract the input message associated to a path
void back_Trellis(connections_node_t Trellis_links[NumStates], decoder_node_t ***Trellis_nodes,int i_state, int i_list, uint8_t (*codeword)[NumSections*NumInputs]){

    int16_t i_section, i_input, i_next_state, initial_state, previous_state, previous_state_bool, actual_state, previous_list_index, i_bit;

    initial_state = i_state;
    actual_state = i_state;
    previous_list_index = (*Trellis_nodes)[NumSections][i_state].previous_index[i_list];
    previous_state_bool = (*Trellis_nodes)[NumSections][i_state].previous_node[i_list];
    previous_state = Trellis_links[i_state].previous_states[previous_state_bool];

    i_bit=NumSections*NumInputs - 1;

    for (i_section = NumSections-1; i_section >= 0 ; i_section--){
        (*codeword)[i_bit]=Trellis_links[actual_state].previous_states_inputs[previous_state_bool][0];
        i_bit--;
        actual_state=previous_state;
        previous_state_bool = (*Trellis_nodes)[i_section][actual_state].previous_node[previous_list_index];
        previous_list_index = (*Trellis_nodes)[i_section][actual_state].previous_index[previous_list_index];
        previous_state = Trellis_links[actual_state].previous_states[previous_state_bool];
    }
    
}

// go backward in the trellis and extract the input message associated to a path
void back_Trellis_small(connections_node_t Trellis_links[NumStates], small_decoder_node_t (*Trellis_nodes)[NumSections+1][NumStates],int i_state, int i_list, uint8_t (*codeword)[NumSections*NumInputs]){

    int16_t i_section, i_input, i_next_state, initial_state, previous_state, previous_state_bool, actual_state, previous_list_index, i_bit;

    initial_state = i_state;
    actual_state = i_state;
    previous_list_index = (*Trellis_nodes)[NumSections][i_state].previous_index[i_list];
    previous_state_bool = (*Trellis_nodes)[NumSections][i_state].previous_node[i_list];
    previous_state = Trellis_links[i_state].previous_states[previous_state_bool];

    i_bit=NumSections*NumInputs - 1;

    for (i_section = NumSections-1; i_section >= 0 ; i_section--){
        (*codeword)[i_bit]=Trellis_links[actual_state].previous_states_inputs[previous_state_bool][0];
        i_bit--;
        actual_state=previous_state;
        previous_state_bool = (*Trellis_nodes)[i_section][actual_state].previous_node[previous_list_index];
        previous_list_index = (*Trellis_nodes)[i_section][actual_state].previous_index[previous_list_index];
        previous_state = Trellis_links[actual_state].previous_states[previous_state_bool];
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


void free2D_Trellis(decoder_node_t **Trellis_nodes){
    int i_row;
    
    for(i_row=0; i_row<=NumSections; i_row++) free(Trellis_nodes[i_row]);
    free(Trellis_nodes);
}






