function [all_TB_messages,d_min,WE]=extract_TB_paths(v,gen_CC,K,d_max)
% This function extracts all tail-biting terminated paths corresponding to
% undetectable errors and whose Hamming distance is smaller than d_max.
% (n,k,v) convolutional code, where 'k' is the number of inputs of the
% convolutional encoder, 'n' the number of outputs, 'v' the memory of the
% encoder.
%
% INPUT(S):
% v : the memory of the convolutional code used (e.g. v=5 if k=1, or v=[3, 7] if k=2, ...)
% gen_CC : generator transfer matrix of the convolutional code in octave
% K : total number of sections of the trellis
% d_max : maximum Hamming distance of the extracted paths
%
% OUTPUT(S):
% all_TB_messages : matrix where each row is the corresponding message of a TB path. The messages are sorted by the distance of the corresponding trellis path.
% d_min : minimum distance of the CC
% WE : weight enumerator up to d_max of the CC. No WE(d=0).
%

%% Parameters
d_max=single(d_max);
num_states = 2^sum(v); % number of trellis states
k=size(gen_CC,1); % number of input bits of the CC encoder
n=size(gen_CC,2); % number of output bits of the CC encoder
M=single(n*K); % high value used to increase the edge weights and remove them from the trellis

% trellis
Trellis=poly2trellis(v+1,gen_CC); % trellis generation
Trellis_nextstate = single(Trellis.nextStates+1);
Trellis_output=zeros(num_states,size(Trellis.outputs,2),'single');
for i_state = 1:num_states
    for j_state = 1:size(Trellis.outputs,2)
        Trellis_output(i_state,j_state)=single(sum(oct2poly(Trellis.outputs(i_state,j_state))));
    end
end
Trellis_output(1,1)= M; % remove the edge S0 -> S0
next_states=2^k; % number of outgoing edges from each state
% sequences of possible output bits
possible_output_bits=generate_all_codewords(eye(k));

%% Initialize outputs
TB_mex = logical([]); % TB terminated messages corresponding to undetectable error events
TB_dist = single([]); % Hamming distance of the found TB messages

for i_state = 1:num_states
    temp_paths=single([i_state,i_state,i_state]);% contains all temporary paths found. [Initial state, Previous State, Current State]
    temp_dist=single(0); % contains the distances of the corresponding temporary paths
    temp_mex = logical([]); % contains the messages of the corresponding temporary paths
    
    TB_length=[]; % length in sections of the found TB messages
        
    for i_section = 1:K
        
        if numel(temp_dist) > 0 % if there are no temporary paths, STOP
            
            temp_paths(:,2) = temp_paths(:,3); % update previous state
        
            % copy previous temporary paths
            c_paths = temp_paths;
            c_lengths = temp_dist;
            c_messages = temp_mex;
        
            % build new temporary paths
            temp_paths = [temp_paths(:,1:2),Trellis_nextstate(temp_paths(:,2),1)];
            temp_dist = [temp_dist+Trellis_output(c_paths(:,2),1)];
            temp_mex=[temp_mex,false(size(c_paths,1),k)];
            for i_nextstate = 2:next_states
                temp_paths = [temp_paths;c_paths(:,1:2),Trellis_nextstate(c_paths(:,2),i_nextstate)];
                temp_dist = [temp_dist;c_lengths+Trellis_output(c_paths(:,2),i_nextstate)];
                temp_mex=[temp_mex;c_messages,repmat(possible_output_bits(i_nextstate,:),size(c_paths,1),1)];
            end
            
            % filter out those paths whose distance is larger than d_max
            index_keep = find(temp_dist<=d_max);
            temp_paths = temp_paths(index_keep,:);
            temp_dist=temp_dist(index_keep,:);
            temp_mex = temp_mex(index_keep,:);
            
            if i_state == 1 % extract undetectable ZT single error events
                
                % save those paths which have merged S0
                index_keep = find((temp_paths(:,1)-temp_paths(:,3)==0));
                TB_dist=[TB_dist;temp_dist(index_keep,:)];
                TB_mex=[TB_mex;[temp_mex(index_keep,:),false(numel(index_keep),k*K-size(temp_mex,2))]];
                TB_length=[TB_length;i_section*ones(numel(index_keep),1)];
                index_keep = find((temp_paths(:,1)-temp_paths(:,3)~=0));

                % filter out those paths which have merged S0
                temp_paths = temp_paths(index_keep,:);
                temp_dist=temp_dist(index_keep,:);
                temp_mex = temp_mex(index_keep,:);
            end
        end
    end

    if i_state == 1 
        d_min = min(TB_dist);
        
        % Concatenate single error events
        % check if there are undetectable single errors which can be concatenated
        if d_max>=2*d_min
            index_find = find(TB_dist <= d_max-d_min);
            % copy undetectable single errors
            one_ZT_mex = TB_mex(index_find,:);
            one_ZT_dist = TB_dist(index_find);
            one_ZT_length = TB_length(index_find);
            % create a list of undetectable multiple errors
            multiple_ZT_mex = one_ZT_mex;
            multiple_ZT_dist = one_ZT_dist;
            multiple_ZT_length = one_ZT_length;
            
            stop_while=0;

            while stop_while<1
                new_ZT_mex=logical([]);
                new_ZT_dist=single([]);
                new_ZT_length=[];
                % concatenate all undetectable single errors with all undetectable
                % multiple errors if the new error has a distance <= d_max and
                % length do not exceed the maximum possible length
                for i_mex = 1 : size(one_ZT_length,1)
                    i_length = one_ZT_length(i_mex);
                    dist_concatenation=one_ZT_dist(i_mex,:)+multiple_ZT_dist;
                    length_concatenation = one_ZT_length(i_mex,:)+multiple_ZT_length;
                    possible_concatenations = find((dist_concatenation<= d_max).*(length_concatenation<=K)).';
                    for j_mex = possible_concatenations
                        j_length = multiple_ZT_length(j_mex);
                        new_dist=dist_concatenation(j_mex);

                        % use the fft correlation to do the circular shift of a message 
                        shifts=-k*[i_length:K-j_length].';
                        all_shifts=abs(ifft((fft(multiple_ZT_mex(j_mex,:)).*exp(1i*2*pi*[0:(K*k)-1].*(shifts)/(K*k))).').')>0.5;

                        % save concatenations
                        new_ZT_mex=[new_ZT_mex;xor(one_ZT_mex(i_mex,:),all_shifts)];
                        new_ZT_dist=[new_ZT_dist;repmat(new_dist,numel(shifts),1)];
                        new_ZT_length=[new_ZT_length;-shifts/k+j_length];
                    end
                end
                % update the TB list
                TB_mex = [TB_mex;new_ZT_mex];
                TB_dist = [TB_dist;new_ZT_dist];
                TB_length = [TB_length;new_ZT_length];
                % update the new undetected multiple errors
                multiple_ZT_mex = new_ZT_mex;
                multiple_ZT_dist = new_ZT_dist;
                multiple_ZT_length = new_ZT_length;

                % check if to stop the while loop
                if numel(multiple_ZT_dist)>0 && min(multiple_ZT_dist)<=d_max-d_min
                    stop_while=0;
                else
                    stop_while=1;
                end  
            end
        end
    end
    
    % remove the i_state from the trellis
    Trellis_output(i_state,:)=M;
    Trellis_output(Trellis_nextstate==i_state)=M;
    
    % if state different from S0, check TB condition of the found paths
    if i_state > 1 && numel(temp_dist) > 1
        index_keep = find((temp_paths(:,1)-temp_paths(:,3)==0));
        TB_dist=[TB_dist;temp_dist(index_keep,:)];
        TB_mex=[TB_mex;temp_mex(index_keep,:)];
    end    
end

% circular shift of all the TB messages
all_TB_messages = TB_mex;
for i_section = 1:K/k
    all_TB_messages = [all_TB_messages;circshift(TB_mex,[0,k*i_section])];
end
TB_dist = repmat(TB_dist,size(all_TB_messages,1)/size(TB_mex,1),1);

% remove equal messages
[all_TB_messages,index_keep,~]=unique(all_TB_messages,'rows');
TB_dist=TB_dist(index_keep,:);

% sort the messages by the distance of their corresponding path
[TB_dist,order]=sort(TB_dist);
all_TB_messages = all_TB_messages(order,:);

% compute the weight enumerator up to d_max, withoud the coefficient of d=0
WE = zeros(1,d_max);
for i_d = 1 : d_max
    WE(i_d)=sum(TB_dist==i_d);
end


