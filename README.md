# Serially Concatenated Convolutional Codes (CC) with outer CRC code
This repository contains codes to
1. design
2. decode

serially concatenated convolutional codes (both tail-biting and zero-tail terminated) with outer CRC codes.

## LICENSE
All the material of this repository is provided under license CC BY-NC-ND

## Author
Riccardo Schiavone, riccardo.schiavone@eurecom.fr

## Related material
Code developed at the German Aerospace Center (DLR) for my Master's Thesis on "Channel Coding for Massive IoT Satellite Systems"
Thesis supervisors: Gianluigi LIVA (DLR), Roberto Garello (Politecnico di Torino), David Gesbert (EURECOM)

## Design best degree-_m_ CRC outer code for a specific CC code
The codes in the folder _Design_CRC_ are written in MATLAB R2018b. They are used to design the best degree-_m_ CRC outer code for a specific convolutional code.

**Inputs**:
- `v` : memory/memories of the CC encoder branch/branches. Same notation of the constraint length of `poly2trellis` function, but this is the memory.
- `gen_CC` : generator polynomials in octave of the convolutional code. Same notation of `poly2trellis` function.
- `K` : length of the information sequence in bits
- `m` : degree of the generator polynomial g(X) of the outer CRC code
- `d_max` : design distance parameter. Attention! The execution time of the code grows exponentially in `d_max`.

**Outputs**:
- `d_min` : minimum distance of the concatenated code
- `best_CRC` : it outputs the generator polynomial of the best CRC code in octave (e.g. x^6+x^5+x^2+1 -> [1 100 101] -> 145)
- `A_min` : it outputs the weight enumerator coefficient at `d_min` of the concatenation of the CC with the `best_CRC` found

The code is based on the following papers:
1. Lou, Chung-Yu, Babak Daneshrad, and Richard D. Wesel. "Convolutional-code-specific CRC code design." IEEE Transactions on Communications 63.10 (2015): 3459-3470.
2. Yang, Hengjie, et al. "An Efficient Algorithm for Designing Optimal CRCs for Tail-Biting Convolutional Codes." 2020 IEEE International Symposium on Information Theory (ISIT). IEEE, 2020.

The algorithms are more efficient than the ones proposed in the previous papers, by the use of the following CRC codes properties:
1. single-bit error detection
2. burst error detection
3. odd detection for those generator polynomials multiples of (X+1)

### Zero-Tail terminated Convolutional Codes (ZTCC)
Run
`[d_min,A_min,best_CRC]=find_best_CRC_ZTCC(v,gen_CC,K,m,d_max)`

### Tail-Biting terminated Convolutional Codes (TBCC)
Run
`[d_min,A_min,best_CRC]=find_best_CRC_TBCC(v,gen_CC,K,m,d_max)`

## Decode a Serially Concatenated CC with an Outer CRC Code
To be added
