function [C] = generate_all_codewords(G)
% This function generates all the codebook C with 2^k codewords of an (n,k)
% linear block code with generator matrix G.
% INPUT(S):
% - G, the binary generator matrix with k independent rows
%
% OUTPUT(S):
% - C, the codebook with 2^k codewords. The data format is a logical matrix.
%
% Attention : if the input is an eye(k) matrix, the output is a matrix with
% all 2^k possible binary vectors with k bits
%

G=logical(G);

C = [false(1,size(G,2));G(1,:)];
for i_row_G = 2 : size(G,1)
    C = [C; xor(G(i_row_G,:),C)];
end
