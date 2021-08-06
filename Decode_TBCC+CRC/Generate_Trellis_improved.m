function [Trellis_matrix, Trellis_state, Trellis_input_matrix]=Generate_Trellis_improved(max_poly_length,generator_polynomial)

% Function to generate the Trellis output matrix and the Trellis input
% matrix of a Convolutional Code with a Generator Polynomial
% each row is a starting state, each column is a final state
%
% INPUT(S) :
% - max_poly_length : maximum degree of the generator polynomial + 1
% - generator_polynomial : can be a vector in octal form (e.g. [5 7] -> '[x²+1 x²+x+1]')
%                          or a string with the polynomials (e.g. {'1 + x^2','1 + x + x^2'})
%
% OUTPUT(S) :
% - Trellis_matrix : it has in each row a starting state, while every output bits columns * index new state, the output bits when go to that new state
%                    If NaN means no connection between those 2 states
%                    ( e.g. Trellis_matrix of a rate 1/2 with 1 bit in input
%                    Trellis_matrix = [ +1, +1, NaN, NaN, -1, -1, NaN, NaN;... % (0,0)
% 				        -1, -1, NaN, NaN, +1, +1, NaN, NaN;... % (0,1)
% 				        NaN, NaN, +1, -1, NaN, NaN, -1, +1;... % (1,0)
% 				        NaN, NaN, -1, +1, NaN, NaN, +1, -1];   % (1,1)
% 				          (0,0) ,  (0,1) ,  (1,0) ,  (1,1)
%                     )
% - Trellis_input_matrix : each row is a starting state, each column the final state, each value is the input symbol needed to reach the final state.
%                          If NaN, no possible transition.

Trellis = poly2trellis (max_poly_length, generator_polynomial) ;

num_out_bits = log2(Trellis.numOutputSymbols);

Trellis_matrix = NaN(Trellis.numStates, size(Trellis.outputs,2) * num_out_bits) ;
Trellis_input_matrix = NaN(Trellis.numStates, Trellis.numStates) ;

for i_state = 1 : Trellis.numStates
    j_column = 1;
    
    for j_state = Trellis.nextStates(i_state,:) % ATTENTION j_state starts from 0
        bits = oct2poly(Trellis.outputs(i_state,j_column));
        if prod(bits==0) == 1
            bits = zeros (1, num_out_bits) ;
        elseif numel(bits) < num_out_bits
            bits = [zeros(1, num_out_bits - numel(bits)) , bits] ;
        end
        Trellis_matrix (i_state, (j_column-1)*num_out_bits+1 : j_column*num_out_bits) = (-1).^(bits) ;
        Trellis_input_matrix (i_state, j_state+1) = j_column - 1 ;

        j_column = j_column + 1 ;
        
        
    end
end

Trellis_state=Trellis.nextStates+1;

				
				
				   
