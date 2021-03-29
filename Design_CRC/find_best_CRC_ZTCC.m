function [d_min,A_min,best_CRC]=find_best_CRC_ZTCC(v,gen_CC,K,m,d_max)
% This function finds the best CRC  wit g(X) of degree m to be concatenated
% with the CC of memory v. It extracts all undetectable errors whose
% trellis path is smaller equal than d_max and check which are undetectable
% errors of the concatenation. Larger d_max, larger the execution time.
%
% INPUT(S):
% v : memory of the CC encoder
% gen_CC : generator polynomials in octave of the convolutional code
% K : length of the information sequence in bits
% m : degree of g(X) of the CRC code
% d_max : design distance parameter
%
% OUTPUT(S):
% d_min : minimum distance of the concatenated codes
% best_CRC : it outputs the best CRC in octave (e.g. x^6+x^5+x^2+1 -> [1 100 101] -> 145)
% A_min : it outputs the weight enumerator coefficient at d_min of the concatenation of the CC with the best_CRC found
%
% License : CC BY-NC-ND
% Author : Riccardo Schiavone, riccardo.schiavone@eurecom.fr
% Developed at the German Aerospace Center (DLR)

%% Extract undetectable single errors of the CC
[ZT_messages,d_min,WE_single,num_shifts,WE]=extract_ZT_singlepaths(v,gen_CC,(K+sum(v)+m)/size(gen_CC,1),d_max);
A_min = WE(d_min);

% property : detect all errors with burst <= degree(CRC generator)
% we remove from the WE all errors with length from the first one to the
% last, which is smaller than m, the degree of the generator polynomial of
% our CRC code
temp_survivors = single(ZT_messages);
temp_survivors(temp_survivors==0)=nan;
temp_survivors = temp_survivors.*[1:size(temp_survivors,2)];
burst_len=max(temp_survivors,[],2,'omitnan')-min(temp_survivors,[],2,'omitnan')+1;
not_surv = find(burst_len <= m);
positions_WE=cumsum(WE_single);
for w=d_min:d_max
    detected_elements = not_surv((not_surv>positions_WE(w-1)) & (not_surv<=positions_WE(w)));
    if numel(detected_elements)>0
        WE_single(w) = WE_single(w) - numel(detected_elements);
        WE(w) = WE(w) - sum(num_shifts(detected_elements));
    end
end
d_min=find(WE_single,1,'first');
surv=find(burst_len > m);
ZT_messages=ZT_messages(surv,:);
num_shifts=num_shifts(surv);
num_messages = size(ZT_messages,1);
length_message = size(ZT_messages,2);

% weight matrix (num_shifts in position (i,j) if mex i has trellis paths weight j, 0 otherwise)
w_ij = false(num_messages,d_max);
start = 1;
stop = 1;
for i_d = d_min : d_max
    start = stop;
    stop = start+WE_single(i_d);
    w_ij(start:stop-1,i_d)=true;
end
w_ij = num_shifts.*w_ij;


%% Search
best_CRC = 2^(m)+1;

cum_WE = cumsum(WE_single);

actual_d = d_min-1;

survived=[2^(m)+1:2:2^(m+1)-1];
H_all_CRC=cell(1,survived(end));
survivors=[];

% check WE of the concatenation with the first 10000 mexages
up_to_d = find(cum_WE<10000,1,'last');

if up_to_d < d_max
    % use the messages of the paths with distance = actual_d
    mex_1=ZT_messages(1:cum_WE(up_to_d),:);
    w_ij_1=w_ij(1:cum_WE(up_to_d),:);
    first_WE = sum(w_ij_1,1);

    % check which messages have even numbers of 1s
    indexes_even_mex=find(mod(sum(mex_1,2),2)==0);
    mex_even=mex_1(indexes_even_mex,:);
    w_ij_even = w_ij_1(indexes_even_mex,:);
    first_WE_even=sum(w_ij_even,1);
else
    mex_1=ZT_messages;
    first_WE = WE;
    w_ij_1=w_ij;
    % check which messages have even numbers of 1s
    indexes_even_mex=find(mod(sum(mex_1,2),2)==0);
    mex_even=mex_1(indexes_even_mex,:);
    w_ij_even = w_ij_1(indexes_even_mex,:);
    first_WE_even=sum(w_ij_even,1);
end

for i_CRC = survived
    gen_CRC = str2num(dec2bin(i_CRC,m+1).').';
    
    [~,rem] = gfdeconv(gen_CRC,[1 1],2);
    
    % check if g(X) of the CRC is multiple of x+1
    % it can detect all messages with odd number of ones
    % if there are no messages with even num of 1s go to next CRC
    if sum(first_WE_even)==0 && sum(rem)==0
        survivors=[survivors,i_CRC];
    else
        G=zeros(length_message-m,length_message);
        G(1,:)=[gen_CRC,zeros(1,length_message-numel(gen_CRC))];
        for i_row = 2:length_message-m
            G(i_row,:)=circshift(G(1,:),i_row-1);
        end
        for i_row = 2:length_message-m
            for j_row = 1:i_row-1
                if G(j_row,i_row)==1
                    G(j_row,:)=mod(G(i_row,:)+G(j_row,:),2);
                end
            end
        end
        H_CRC=[G(:,length_message-m+1:end);eye(m)];
        H_all_CRC{i_CRC}=logical(H_CRC);
        % if g(X) no multiple of x+1, check all messages whose path
        % has distance = actual_d, otherwise only those messages
        % with even number of 1s
        if sum(rem)>0
            syndromes = (sum(mod(mex_1*H_CRC,2),2)>0).';
            WE_actual_d = first_WE - sum(w_ij_1(syndromes,:),1);
        else
            syndromes = (sum(mod(mex_even*H_CRC,2),2)>0).';
            WE_actual_d = first_WE_even - sum(w_ij_even(syndromes,:),1);
        end
        
        
        % check if the CRC code can detect all paths with weight = actual_d
        if sum(WE_actual_d)==0
            survivors=[survivors,i_CRC];
        else
            d_min_CRC = find(WE_actual_d,1,'first');
        end
        
        % if there is no CRC code capable of detect all paths with
        % Hamming weight = actual_d, then save the best
        if numel(survivors)==0
            if d_min_CRC > d_min
                d_min = d_min_CRC;
                best_CRC = i_CRC;
                A_min = WE_actual_d(d_min_CRC);
            elseif d_min_CRC == d_min
                if WE_actual_d(d_min_CRC) < A_min
                    best_CRC = i_CRC;
                    A_min = WE_actual_d(d_min_CRC);
                end
            end
        end
    end
end

actual_d = up_to_d;

% check WE of the remaining messages
% if there are any survived CRC candidates
while numel(survivors)>0
    actual_d=actual_d+1;
    if actual_d > d_max
        error("d_max insufficient!");
    end
    
    if start==1
        start=0;
    else
        survived=survivors;
        survivors=[];
    end
    
    % use the messages of the paths with distance = actual_d
    mex_1=ZT_messages(cum_WE(actual_d-1)+1:cum_WE(actual_d),:);
    first_WE = WE(actual_d);
    num_shifts1=num_shifts(cum_WE(actual_d-1)+1:cum_WE(actual_d));
    % check which messages have even numbers of 1s
    indexes_even_mex=find(mod(sum(mex_1,2),2)==0);
    mex_even=mex_1(indexes_even_mex,:);
    num_shifts_even = num_shifts1(indexes_even_mex,:);
    first_WE_even=sum(num_shifts_even);
    
    % if no messages of paths with Hamming weight = actual_d, go to next weight
    if numel(mex_1)==0
        survivors=survived;
    else
        
        % check how many messages whose corresponding codewords have
        % Hamming weight == actual_d can be detected by the survived CRC codes
        for i_CRC = survived
            gen_CRC = str2num(dec2bin(i_CRC,m+1).').';
            
            [~,rem] = gfdeconv(gen_CRC,[1 1],2);
            
            % check if g(X) of the CRC is multiple of x+1
            % it can detect all messages with odd number of ones
            % if there are no messages with even num of 1s go to next CRC
            if first_WE_even==0 && sum(rem)==0
                survivors=[survivors,i_CRC];
            else
                
                H_CRC=H_all_CRC{i_CRC};
                
                % if g(X) no multiple of x+1, check all messages whose path
                % has distance = actual_d, otherwise only those messages
                % with even number of 1s
                if sum(rem)>0
                    syndromes = (sum(mod(mex_1*H_CRC,2),2)>0).';
                    WE_actual_d = first_WE - syndromes*num_shifts1;
                else
                    syndromes = (sum(mod(mex_even*H_CRC,2),2)>0).';
                    WE_actual_d = first_WE_even - syndromes*num_shifts_even;
                end
                
                % check if the CRC code can detect all paths with weight = actual_d
                if WE_actual_d==0
                    survivors=[survivors,i_CRC];
                else
                    d_min_CRC = actual_d;
                end
                
                % if there is no CRC code capable of detect all paths with
                % Hamming weight = actual_d, then save the best
                if numel(survivors)==0
                    if d_min_CRC > d_min
                        d_min = d_min_CRC;
                        best_CRC = i_CRC;
                        A_min = WE_actual_d;
                    elseif d_min_CRC == d_min
                        if WE_actual_d < A_min
                            best_CRC = i_CRC;
                            A_min = WE_actual_d;
                        end
                    end
                end
            end
        end
    end
end
g_CRC = dec2base(bin2dec(num2str(fliplr(str2num(dec2bin(best_CRC,m+1).').'))),8);
fprintf('\nConcatenation with the best CRC code with degree m=%d has\nd_min=%d, A_min=%d and g_CRC(X)=%s',m,d_min,A_min,g_CRC)
fprintf(" (e.g. g(X) = x^4+x+1 = [1 0 0 1 1] = 23)\n")

