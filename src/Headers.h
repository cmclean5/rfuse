#if !defined(HEADERS_INCLUDED)
#define HEADERS_INCLUDED

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(bigstatsr, rmio)]]
// [[Rcpp::depends(RcppArmadillo)]]

//C and C++
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <iostream>
#include <string>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <cstdint>     //defines uint_32, uint64 etc...
#include <cstddef>     //defines std::size_t (unsigned) and std::ssize_t (signed)
#include <stdlib.h>    // to use malloc and free
//#include <x86intrin.h> // for SSE intinstics, gcc ver. >=8 (https
#include <mm_malloc.h>

//OpenMP header
#include <omp.h>

//Rcpp
//#include <Rcpp.h>

//Armadillo
#include <RcppArmadillo.h>
//#include <armadillo>

//bigstatsr
#include <bigstatsr/BMAcc.h>


//Name space
using namespace arma;
using namespace Rcpp;
using namespace std;

//*** CONSTANTS ***

//define constants
#define CSIZE  64 //the cache line size (NOTE: can find the cache line size using: getconf LEVEL1_DCACHE_LINESIZE  , for 'big' this is 64. 

#define CHUNK  100    //for openmp scheduler's
#define BLOCK_SIZE 16
#define SMALL  1.0e-100
#define ROUND_UP(x, s) (((x)+((s)-1)) & -(s))

#endif
