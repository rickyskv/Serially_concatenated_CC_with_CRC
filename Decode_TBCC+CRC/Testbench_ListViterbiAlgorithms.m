%% Testbench to measure the performances of the Parallel-List and Serial-List Viterbi Algorithm
clear variables
close all

%% Parameters
EbN0=[1:0.5:3]; % in dB
decoder='SLVA'; % 'PLVA' or 'SLVA'
k=64; % length of the info message
num_errors = 500; % number of errors counted for each simulated EbN0 value
num_WAVA_iterations = 1; % number of WAVA iterations (1 is the simple Viterbi Algorithm)
num_processes = 8; % to speed-up simulations. It depends by your machine
g_CC = [133,171]; % generator polynomials of the convolutional encoder in octave
g_CRC = 177; % generator polynomial of the CRC code in octave (g(x) = 1+x^2+x^4+x^5 -> g = [1 0 1 0 1 1]_2 -> [53]_8 )
max_list_size = [1];%,4,16,32];%[1,1e2,1e4]; % vector of the various maximum list sizes to be tested (e.g. [1,1e2,1e4] for SLVA, [1,8,64] for PLVA)
compiler='gcc'; % C compiler (e.g. gcc, icc, ...). Note that only linux compiler are compatible
optimizer='-O3'; % optimization flag of the compiler (e.g. -O3, -O0, -O1, -Ofast, ...)
axis_values=[EbN0(1),EbN0(end)+1e-10,1e-4,1]; % axis values for plots

v=numel(oct2poly(g_CC(1)))-1; % memory of the convolutional encoder
m=numel(oct2poly(g_CRC))-1; % degree of the polynomial of the outer CRC code
R = 1/size(g_CC,2); % rate of the convolutional encoder

% the desired puncturing vector (1 if element not punctured, 0 if punctured)
punct_v = ones(1,(k+m)/R);
punct_v(11:11:end) = 0;

% strings
if not(isfolder('Results'))
    mkdir('Results');
end
if not(isfolder('Inputs'))
    mkdir('Inputs');
end
out_file_header = strcat('./Results/Test_1_',decoder,'_'); % header of the output file with the results
title_string=''; % title of the plot
New_T_file = 'trellis.h';
puncturing_vector_file = './Inputs/puncturing.txt'; % where to save the puncturing vector

%% Useful variables and structures
colors='gbrcmgbrcm';

if strcmp(decoder,'PLVA')
    C_file='New_PLVA_predecessor';
    decoder_string = 'PARALLEL-LIST VITERBI ALGORITHM';
elseif strcmp(decoder,'SLVA')
    C_file='New_SLVA_predecessor';
    decoder_string = 'SERIAL-LIST VITERBI ALGORITHM';
end

if numel(EbN0)==1
    stepsize=1e-4;
else
    stepsize = EbN0(2)-EbN0(1);
end


num_outs=1/R; % output bits of the convolutional encoder for each input bits
n=k/R; % nominal number of output bits from convolutional encoder without outer code

% Generate trellis structure for (1/R,1,v) TBCC encoders
[Trellis, Trellis_nextstate, Trellis_input_symb]=Generate_Trellis_improved(v+1,g_CC);

Print_puncturing_vector(punct_v, puncturing_vector_file);
if sum(1-punct_v) > 0
    name_out = int2str(sum(1-punct_v));
else
    name_out='no';
end

g_TBCC=zeros(numel(g_CC),numel(oct2poly(g_CC(1))));
for i=1:numel(g_CC)
    g_TBCC(i,:) = oct2poly(g_CC(i));
end

if m > 0
    name_CRC = 'yes';
else
    name_CRC = 'no';
end


%% Simulations
fprintf("\n------ %s ------\n",decoder_string)
fprintf("\n(%d,%d) TBCC+CRC code\n",(k+m)/R-sum(1-punct_v),k)
fprintf("TBCC: memory %d, generators [%d,%d]_8\nCRC: degree %d, polynomial [%d]_8\n", v, g_CC(1), g_CC(2), m, g_CRC)
fprintf("\n%d Wrap-Around Viterbi Algorithm iterations\n",num_WAVA_iterations-1)

figure('units','normalized','outerposition',[0 0 1 1])
hold on

for list_size = max_list_size
	fprintf("\n--- Max list size L=%d\n\n",list_size)

	OUT_file = strcat(out_file_header,'_memory',int2str(sum(v)),name_CRC,'_CRC_',name_out,'_puncturing_L',int2str(list_size),'_numWAVA',int2str(num_WAVA_iterations),'.txt');
	CER = zeros(1,numel(EbN0)); % codeword error rate

	%% Print Trellis and compile the C file
	if strcmp(decoder,'PLVA')
	    New_Print_Trellis(Trellis, Trellis_nextstate, Trellis_input_symb, New_T_file, k, m, max(2,list_size), g_TBCC, oct2poly(g_CRC));
	elseif strcmp(decoder,'SLVA')
	    New_Print_Trellis(Trellis, Trellis_nextstate, Trellis_input_symb, New_T_file, k, m, 2, g_TBCC, oct2poly(g_CRC));
	end
	system('rm -f ./*.o');
	system(strjoin({'gcc',strcat(New_T_file(1:end-1),'c'),strcat(C_file,'.c'),'-o',strcat(C_file,'.o'),'-lm',optimizer,'-w'}));


    
    %% Simulation
    start=tic;
    
    % string of the command
     c_command=strjoin({strcat('./',C_file,'.o'), OUT_file, int2str(num_WAVA_iterations), int2str(list_size), int2str(num_errors), puncturing_vector_file, int2str(num_processes), num2str(EbN0(1),'%.1f'), num2str(EbN0(end),'%.1f'), num2str(stepsize,'%.5f')});

    % run the simulation in C
    system(c_command);
    
    total_time=toc(start);
    
    fprintf("\nExecution time: %.2f seconds",total_time)

    %% Read results
    CERs=readtable(OUT_file,'ReadVariableNames',false);
    CER=cell2mat(CERs.Var7);
    CER=str2num(CER(:,1:end-1));
    total_simulated_messages = floor(sum(ceil(num_errors/num_processes)*num_processes./CER));
    fprintf("\nNumber of simulated messages to find %d errors:\n%d messages (%d messages/second)\n",ceil(num_errors/num_processes)*num_processes,total_simulated_messages,floor(total_simulated_messages/total_time))
    
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
semilogy([1:0.5:3],[1e-1,4e-2,8.5e-3,1.1e-3,1e-4],'--k','DisplayName','RCU bound (128,64)')

title(title_string)
xlabel('E_b / N_0 [dB]'),ylabel('CER'), grid on
axis(axis_values)
set(gca, 'YScale', 'log')
legend('show','NumColumns',4,'Location','southoutside','Orientation','horizontal')

print(strcat(out_file_header,'v',int2str(v),'_m',int2str(m),'_(',int2str((k+m)/R-sum(1-punct_v)),',',int2str(k),')_I',int2str(num_WAVA_iterations),'_L',sprintf('_%d',max_list_size)),'-depsc')

