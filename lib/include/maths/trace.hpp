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

#include "maths/arma.hpp"
#include "utils/std.hpp"
#include <stdexcept>

namespace maths{
using namespace maths::armadillo;
using namespace utils::stds;
using std::invalid_argument;

/*
 * The trace function performs a special kind of summation: it sums sub-matrix
 * blocks of a matrix along the diagonal. Works for any armadillo matrix 
 * data structures.
 *
 * -----------------------------------------------------------------------------
 * A --------> Square matrix over which the trace is performed.
 * N --------> Size of the sub-matrix.
 * vertices -> Index of diagonal matrices that are summed over.
 * -----------------------------------------------------------------------------
 */

template<class T>
T trace(const T &A, uint N = 1, const ucol *vertices = 0){
    int m = A.n_rows;
    int n = A.n_cols;
    
    if (m != n){
        stringstream err;
        err << "In trace(A, N): A(" << m << "x" << n << ")  is not square matrix.";
        throw invalid_argument(err.str());
    }    
    if (m%N != 0){
        stringstream err;
        err << "In trace(A, N): number of rows/cols of A(" << m << "x" << n << ") has to be multiple of N = " << N << ".";
        throw invalid_argument(err.str());
    } 
    if(vertices != 0 && (max(*vertices)+1)*N-1 > m){
        stringstream err;
        err << "In trace(A, N, vertices): one of the vertices exceeds size of of A(" << m << "x" << n << "). max(vertices)*(N+1)-1 = " <<  max(*vertices)*(N+1)-1 << ".";
        throw invalid_argument(err.str());        
    }
    
    if (m == N){
        return A;
    }

    T tr = zeros<T>(N, N);
    
    // optimized for speed. requires more code.
    if(vertices == 0){// trace over all elements
        for(uint i = 0; i < m; i += N){
            uint j = i+N-1;
            tr = tr + A(span(i,j),span(i,j));
        }
    }else{// trace over selective elements
        const ucol& vert = *vertices;
        for(uint i = 0; i < vert.n_elem; ++i){
            uint beg = vert[i]*N;
            uint end = beg + N - 1;
            tr = tr + A(span(beg,end),span(beg,end));
        }        
    }
    return tr;
}

}
#endif	/* TRACE_HPP */

