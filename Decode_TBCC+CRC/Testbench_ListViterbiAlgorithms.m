%% Testbench to measure the performances of the Parallel-List Viterbi Algorithm
clear variables
close all

%% Parameters
EbN0=[1.0:0.5:3.0]; % in dB
decoder='SLVA'; % 'PLVA' or 'SLVA'
k=64; % length of the info message
R=1/2; % rate of the Convolutional Code
num_errors = 100; % number of errors counted for each simulated EbN0 value
num_WAVA_iterations = 1; % number of WAVA iterations (1 is Viterbi)
num_processes = 20; % to speed-up simulations. It depends by your machine
v = 6; % memory of the convolutional encoder
g_CC = [133,171]; % generator polynomials of the convolutional encoder in octave
m = 6; % degree of the polynomial of the outer CRC code
g_CRC = 101; % generator polynomial of the CRC code in octave
max_list_size = [1,100,1e4]; % vactor of the various list sizes to test
punctures = m; % number of punctured input bits
compiler='gcc'; % C compiler (e.g. gcc, icc, ...). Note that only linux compiler are compatible
optimizer='-O3'; % optimization flag of the compiler (e.g. -O3, -O0, -O1, -Ofast, ...)
axis_values=[EbN0(1),EbN0(end),1e-4,1]; % plot axis values

% strings
out_file_header = strcat('Test_1_',decoder,'_'); % header of the output file with the results
title_string=''; % title of the plot
H_file = 'H_CRC.txt'; % where to save the parity check matrix
T_file = 'trellis.txt'; % where to save the trellis structure

%% Useful variables and structures
colors='gbrcmgbrcm';

if strcmp(decoder,'PLVA')
    C_file='ParallelListViterbiAlgorithm';
    decoder_string = 'PARALLEL-LIST VITERBI ALGORITHM';
elseif strcmp(decoder,'SLVA')
    C_file='SerialListViterbiAlgorithm';
    decoder_string = 'SERIAL-LIST VITERBI ALGORITHM';
end

num_outs=1/R; % output bits of the convolutional encoder for each input bits
n=k/R; % nominal number of output bits from convolutional encoder without outer code

% Generate trellis structure for (1/R,1,v) TBCC encoders
[Trellis, Trellis_nextstate, Trellis_input_symb]=Generate_Trellis_improved(v+1,g_CC);

% Print trellis in a file
Print_Trellis(Trellis, Trellis_nextstate, Trellis_input_symb, T_file);

% Generate the generator matrix of the TBCC
u = [1, zeros(1,k+m-1)];
G_CC = zeros(k+m,(k+m)/R);
T=poly2trellis(v+1,g_CC);
G_CC(1,:) = convenc(u,T);
for i=2:k+m
    G_CC(i,:)=circshift(G_CC(1,:),(1/R)*(i-1));
end

% generate the parity check matrix of the CRC code for syndrome check
% if memory m=0, the parity check matrix is a allzero vector, which means
% that all syndrome checks are ok
gen_CRC=oct2poly(g_CRC);
if m>0
    G=zeros(k,k+m);
    % generate G
    G(1,:)=[gen_CRC,zeros(1,k-1)];
    for i_row = 2:k
        G(i_row,:)=circshift(G(1,:),i_row-1);
    end
    % make it systematic
    for i_row = 2:k
        for j_row = 1:i_row-1
            if G(j_row,i_row)==1
                G(j_row,:)=mod(G(i_row,:)+G(j_row,:),2);
            end
        end
    end
    % systematic H
    H_CRC=[G(:,k+1:end);eye(m)].';
else
    H_CRC=zeros(1,k);
end

if punctures > 0
    name_out = int2str(punctures);
else
    name_out='no';
end

if m > 0
    name_CRC = 'yes';
else
    name_CRC = 'no';
end

Print_CRC(H_CRC,H_file);

G=[eye(k)];
if m> 0
    G_CRC = generator_CRC(size(G,2),m,gen_CRC);
    G = mod(G * G_CRC,2);
end

G = mod(G * G_CC,2);

G_file=strcat('G_v',int2str(v),'_m',int2str(m),'.txt');

Print_G(G,G_file);

%% Simulations
fprintf("\n------ %s ------\n",decoder_string)
fprintf("\n(%d,%d) TBCC+CRC code\n",(k+m-punctures)/R,k)
fprintf("TBCC: memory %d, generators [%d,%d]_8\nCRC: degree %d, polynomial [%d]_8\n", v, g_CC(1), g_CC(2), m, g_CRC)

figure('units','normalized','outerposition',[0 0 1 1])
hold on

for list_size = max_list_size
    fprintf("\n--- List size L=%d\n\n",list_size)
    
    OUT_file = strcat(out_file_header,'_memory',int2str(sum(v)),name_CRC,'_CRC_',name_out,'_puncturing_L',int2str(list_size),'_numWAVA',int2str(num_WAVA_iterations),'.txt');
    CER = zeros(1,numel(EbN0)); % codeword error rate
    
    %% Compile the C file at the first simulated list size and remove old compiled files
    if list_size == max_list_size(1)
        system('rm -f ./*.o');
        system(strjoin({'gcc',strcat('./',C_file,'.c'),'-o',strcat(C_file,'.o'),'-lm',optimizer,'-w'}));
    end
    
    %% Simulation
    start=tic;
    % string of the command
    c_command=strjoin({strcat('./',C_file,'.o'),T_file,H_file,OUT_file,int2str(num_WAVA_iterations),int2str(list_size),int2str(num_errors),int2str(k),int2str(k/R),int2str(m),int2str(punctures),int2str(num_processes),G_file,num2str(EbN0(1),'%.1f'),num2str(EbN0(end),'%.1f'),num2str(EbN0(2)-EbN0(1),'%.1f')});

    % run the simulation in C
    system(c_command);
    total_time=toc(start);
    fprintf("\nExecution time: %.2f seconds",total_time)

    %% Read results
    CER=readtable(OUT_file,'ReadVariableNames',false);
    CER=CER.Var1;
    total_simulated_messages = floor(sum(ceil(num_errors/num_processes)*num_processes./CER));
    fprintf("\nNumber of simulated messages to find %d errors:\n%d messages (%d messages/second)\n",num_errors,total_simulated_messages,floor(total_simulated_messages/total_time))
    
    %% Plot results
    legend_label=strcat('\nu=',int2str(v),'+CRC-',int2str(m),', I=',int2str(num_WAVA_iterations),', L=',int2str(list_size));
    i_list = find(max_list_size == list_size);
    i_v = find(v == [3,6,8]);
    if num_WAVA_iterations==1
        semilogy(EbN0,CER,'-.o','DisplayName',legend_label,'Color',colors(i_list));
    elseif num_WAVA_iterations==2
        semilogy(EbN0,CER,'-d','DisplayName',legend_label,'Color',colors(i_list));
    elseif num_WAVA_iterations==3
        semilogy(EbN0,CER,'-s','DisplayName',legend_label,'Color',colors(i_list));
    else
        semilogy(EbN0,CER,'-|','DisplayName',legend_label,'Color',colors(i_list));
    end

end

%% Plot known results for comparison and set figure properties
semilogy([1:0.5:3],[1e-1,4e-2,8.5e-3,1.1e-3,1e-4],'--k','DisplayName','RCU bound')

title(title_string)
xlabel('E_b / N_0 [dB]'),ylabel('CER'), grid on
axis(axis_values)
set(gca, 'YScale', 'log')
legend('show','NumColumns',4,'Location','southoutside','Orientation','horizontal')

print(strcat('v',int2str(v),'_m',int2str(m),'_',decoder),'-depsc')

