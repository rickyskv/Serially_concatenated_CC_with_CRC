# Serially Concatenated Convolutional Codes (CC) with outer CRC code
This repository contains codes to
1. design
2. decode

serially concatenated convolutional codes (both tail-biting and zero-tail terminated) with outer CRC codes.

## LICENSE
All the material of this repository is provided under license CC BY-NC-ND

## Author
Riccardo Schiavone, riccardo.schiavone@polito.it, rs.riccardoschiavone@gmail.com

## Related material
Code developed at the German Aerospace Center (DLR) for my Master's Thesis on "Channel Coding for Massive IoT Satellite Systems"
Thesis supervisors: Gianluigi LIVA (DLR), Roberto Garello (Politecnico di Torino), David Gesbert (EURECOM)

## 1) Decode serially concatenated tail-biting convolutional codes with outer CRC code
*Note: the code is compatible with linux machines. Multi-process configuration is possible, simply changing the parameter in the testbench file. If you are not using a machine with an Intel processor, remove the libraries in the C files used for vector computations.*

Files:
- `Testbench_ListViterbiAlgorithms.m` generates the files (trellis, generator matrix, parity check matrix) used by the simulators implemented in the C files. Modify the parameters section to test various codes.
- `ParallelListViterbiAlgorithms.c` implements the TBCC+CRC encoder, adds AWGN impairments, punctures the code, and decodes the received vectors with the given parameters.
- `SerialListViterbiAlgorithms.c` implements the TBCC+CRC encoder, adds AWGN impairments, punctures the code, and decodes the received vectors with the given parameters.
- `Print_*.m` prints the matrices/structures used by the simulators implemented in the C files.
- `generator_*.m` and `Generate_*.m` generate some useful matrices needed in the code.

In case of any errors or problems, please contact me via email (riccardo.schiavone@polito.it, rs.riccardoschiavone@gmail.com)

## 2) Design best degree-_m_ CRC outer code for a specific CC code
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

If `d_max` is not sufficient to find the best degree-_m_ CRC code, the code will display to your screen an error message. In that case run again your algorithm with a larger `d_max`.

The code is based on the following papers:
1. Lou, Chung-Yu, Babak Daneshrad, and Richard D. Wesel. "Convolutional-code-specific CRC code design." IEEE Transactions on Communications 63.10 (2015): 3459-3470.
2. Yang, Hengjie, et al. "An Efficient Algorithm for Designing Optimal CRCs for Tail-Biting Convolutional Codes." 2020 IEEE International Symposium on Information Theory (ISIT). IEEE, 2020.

The algorithms are more efficient than the ones proposed in the previous papers in the number of checks, by the use of the following CRC codes properties:
1. single-bit error detection
2. burst error detection
3. odd detection for those generator polynomials multiples of (X+1)

### Zero-Tail terminated Convolutional Codes (ZTCC)
The algorithm will consider a _(K+m+v,K)_ serially concatenated code.

Run: 
`[d_min,A_min,best_CRC]=find_best_CRC_ZTCC(v,gen_CC,K,m,d_max);`

e.g. to find the best degree-6 CRC code for the famous constraint length 7 ZTCC with generator polynomial [133,171], when `K=64`, run: 

`[d_min,A_min,best_CRC]=find_best_CRC_ZTCC(7-1,[133,171],64,6,14);`

### Tail-Biting terminated Convolutional Codes (TBCC)
The algorithm will consider a _(K+m,K)_ serially concatenated code.

Run: 
`[d_min,A_min,best_CRC]=find_best_CRC_TBCC(v,gen_CC,K,m,d_max);`

e.g. to find the best degree-6 CRC code for the famous constraint length 7 ZTCC with generator polynomial [133,171], when `K=64`, run: 

`[d_min,A_min,best_CRC]=find_best_CRC_ZTCC(7-1,[133,171],64,6,14);`

## Decode a Serially Concatenated CC with an Outer CRC Code
To be added
