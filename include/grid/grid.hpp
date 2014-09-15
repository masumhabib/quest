/* 
 * File:   grid.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 26, 2014, 10:08 AM
 */

#ifndef GRID_HPP
#define	GRID_HPP

#include "utils/Printable.hpp"
#include "maths/svec.h"
#include <armadillo>
#include <iostream>
#include <limits>

namespace utils{
using namespace maths::armadillo;
using namespace maths::spvec;
using std::stringstream;
using std::numeric_limits;

/*
 * Generates linear grid
 */
template<class T>
Col<T> linspace(T start, T end, T delta);
template<class T>
void linspace(Col<T>& col, T start, T end, T delta);

/*
 * Mesh grid: two dimensional, linear matrix like grid
 * 
 *  X = | 1 2 3 |      Y = | 6 6 6 |
 *      | 1 2 3 |          | 5 5 5 |
 *  iy >| 1 2 3 |      iy >| 4 4 4 |
 *        ^                  ^
 *       ix                  ix
 */
template<class T>
void meshgrid(Mat<T> &X, Mat<T> &Y, const Col<T> &xi, const Col<T> &yi);
template<class T>
void meshgrid(Col<T> &X, Col<T> &Y, const Col<T> &xi, const Col<T> &yi);


/* Definitions */
template<class T>
Col<T> linspace(T start, T end, T delta){
    long N;
    T correctEnd;
    if (delta <= numeric_limits<T>::min()){
        N = 1;
    }else{
        N = floor(abs((end - start)/delta))+1;
        correctEnd = start + (N - 1)*delta;
    }
    
    return maths::armadillo::linspace<Col<T> >(start, end, N);
}

template<class T>
void linspace(Col<T>& col, T start, T end, T delta){
    col = linspace<T>(start, end, delta);
}

template<class T>
void meshgrid(Mat<T> &X, Mat<T> &Y, T startx, T endx, T delx, T starty, T endy, T dely){
    Col<T> x = linspace<T>(startx, endx, delx);
    Col<T> y = linspace<T>(starty, endy, dely);

    meshgrid<T>(X, Y, x, y);         
}

template<class T>
void meshgrid(Col<T> &X, Col<T> &Y, T startx, T endx, T delx, T starty, T endy, T dely){
    Col<T> x = linspace<T>(startx, endx, delx);
    Col<T> y = linspace<T>(starty, endy, dely);

    meshgrid<T>(X, Y, x, y);         
}

template<class T>
void meshgrid(Mat<T> &X, Mat<T> &Y, const Col<T> &xi, const Col<T> &yi){
    int NX = xi.n_elem;
    int NY = yi.n_elem;
    
    X.set_size(NY, NX);
    Y.set_size(NY, NX);

    for(int iy = 0; iy < NY; ++iy){
        for(int ix = 0; ix < NX; ++ix){
            X(iy,ix) = xi(ix);
            Y(iy,ix) = yi(iy);
        }
    }
}

template<class T>
void meshgrid(Col<T> &X, Col<T> &Y, const Col<T> &xi, const Col<T> &yi){
    
    Mat<T> matx;
    Mat<T> maty;
    
    meshgrid(matx, maty, xi, yi);
    
    X = vectorise(matx);
    Y = vectorise(maty);
}

}
#endif	/* GRID_HPP */

