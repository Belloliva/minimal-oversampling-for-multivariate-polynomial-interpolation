C code for producing the numerical experiemts of the paper 

**Chkifa, M.A, Dolbeault, M. 
Randomized least-squares with minimal oversampling and interpolation in general spaces.**

To appear in SIAM Journal on Numerical Analysis (SINUM)

The code requires to have gnu scientific library (gsl) installed. 

To compile and run (lunix/macos): 

gcc -Wall -w -I/usr/local/include -c *.c

gcc -L/usr/local/lib *.o -lgsl -lgslcblas 

./a.out


Structure of the code: the files are part of a multivariate polynomial approximation (mpa) library. 
The files implements the concepts of mutliindex, multivariate polynomial, lower set, etc 

- mpa_defines.h : simple macros
- mpa_combinatorics.h /.c : simple combinatorics functions 
- mpa_stdlib.h /.c : few functions for inserting in ordered array   
- mpa_multi_index.h /.c: multiindex notation 
- mpa_monotone_graph.h /.c: datastructure for lower set(also called downward closed) 
- mpa_polynomial.h /.c: usual polynomials
- mpa_function.h : 
- mpa_sampling.h /.c: sampling for weighted least squares, the algorithm of the paper included  
- mpa_linalg.h /.c: few simple linear algebra functions 
- Least_squares_solver.h/.c: solver for computing the least square projection 
- main.c


