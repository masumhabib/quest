/* 
 * File:   trace.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 23, 2014, 4:31 PM
 * 
 * Implementation of generalized trace function.
 * 
 */

#ifndef TRACE_HPP
#define	TRACE_HPP

#include "arma.hpp"
#include <stdexcept>

using namespace maths::armadillo;
using std::invalid_argument;

/*
 * The trace function performs a special kind of summation: it sums sub-matrix
 * blocks of a matrix along the diagonal. Works for any armadillo matrix 
 * data structures.
 *
 * -----------------------------------------------------------------------------
 * A --------> Square matrix over which the trace is performed.
 * N --------> Size of the sub-matrix.
 * -----------------------------------------------------------------------------
 */

template<class T>
T trace(const T &A, uint N = 1){
    int m = A.n_rows;
    int n = A.n_cols;
    
    if (m != n){
        throw invalid_argument("In trace(A, N), A is not square matrix.");
    }    
    if (m%N != 0){
        throw invalid_argument("In trace(A, N), number of elements of A has to be multiple of N");
    } 
    
    if (m == N){
        return A;
    }

    T tr = zeros<T>(N, N);
    for(int i = 0; i < m; i += N){
        tr = tr + A(span(i,i+N-1),span(i,i+N-1));
    }
    
    return tr;
}

#endif	/* TRACE_HPP */

