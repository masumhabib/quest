/* 
 * File:   grid.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 26, 2014, 10:08 AM
 */

#ifndef GRID_HPP
#define	GRID_HPP

#include "../utils/Printable.hpp"
#include "../maths/svec.h"
#include <armadillo>
#include <iostream>

namespace utils{
using namespace maths::armadillo;
using namespace maths::spvec;
using std::stringstream;

template<class T>
void meshgrid(Mat<T> &X, Mat<T> &Y, const Col<T> &xi, const Col<T> &yi);
template<class T>
void meshgrid(Col<T> &X, Col<T> &Y, const Col<T> &xi, const Col<T> &yi);

/*
 * Vector grid: a one dimensional, linear grid.
 */
template<class T>
struct Grid1D:public Printable{
public:   
protected: 
    // data structure
    T           mMinV;
    T           mMaxV;
    T           mdV;
    int         mNV;
    Col<T>      mV;
    
public:    
    Grid1D(Col<T>& V, string prefix = ""):Printable(" " + prefix), mV(V){
        mTitle = "Grid";

        mNV = V.n_elem;
        mMaxV = arma::max(V);
        mMinV = arma::min(V);

        mdV = mV(1)-mV(0);
    }

    Grid1D(T min, T max, T d, string prefix = ""):Printable(" " + prefix),
        mMinV(min), mMaxV(max), mdV(d)
    {
        mTitle = "Grid";

        if (mdV == 0){
            mNV = 1;
        }else{
            mNV = floor(abs((mMaxV - mMinV)/mdV))+1;
            mMaxV = mMinV + (mNV - 1)*mdV;
        }
        create();
        //mMaxV = mV(mNV-1);        
    }

    Grid1D(T min = 0, T max = 0, int N = 1, string prefix = ""):
        Printable(" " + prefix), mMinV(min), mMaxV(max), mNV(N)
    {
        mTitle = "Grid";

        create();
        if (mNV == 1){
            mdV = 0;
        }else{
            mdV = mV(1)-mV(0);
        }
    }
    
    virtual T operator()(int i){
        return mV(i);
    } 
    
    const Col<T>& operator()(){
        return mV;
    } 
    
    T   min() const { return mMinV; };
    T   max() const { return mMaxV; };
    T   del() const { return mdV; };
    int N()   const { return mNV; };
    
    T   V(int it) { return (*this)(it); }
    
    string toString() const { 
        stringstream ss;
        ss << Printable::toString() << ":" << endl;
        ss << mPrefix << " min = " << min()
                      << ", max = " << max() 
                      << ", N = " << N() << ", del = " << del() << endl;
        ss << mPrefix << mV;

        return ss.str(); 
    }; 
    
protected:
    virtual void create(){
        mV = linspace<Col<T> >(mMinV, mMaxV, mNV);
    }
};

/*
 * Matrix grid: two dimensional, linear matrix like grid
 * 
 *  X = | 1 2 3 |      Y = | 6 6 6 |
 *      | 1 2 3 |          | 5 5 5 |
 *  iy >| 1 2 3 |      iy >| 4 4 4 |
 *        ^                  ^
 *       ix                  ix
 * 
 */

template<class T>
struct Grid2D:public Printable{
public:
protected:
    // data structure
    T           mMinX;
    T           mMaxX;
    T           mdX;
    T           mNX;
    T           mMinY;
    T           mMaxY;
    T           mdY;
    T           mNY;
    
    Mat<T>      mX;
    Mat<T>      mY;
    
public:
    
    Grid2D(T minx, T maxx, T dx, T miny, T maxy, T dy, string prefix = ""):
        Printable(" " + prefix),
        mMinX(minx), mMaxX(maxx), mdX(dx), mMinY(miny), mMaxY(maxy), mdY(dy)
    {
        mTitle = "2D Grid";

        if (mdX == 0){
            mNX = 1;
        }else{
            mNX = abs((mMaxX - mMinX)/mdX)+1;
        }
        if (mdY == 0){
            mNY = 1;
        }else{
            mNY = abs((mMaxY - mMinY)/mdY)+1;
        }
        create();
        mMaxX = mX(mNX-1);
        mMaxY = mY(mNY-1);
    }

    Grid2D(T minx, T maxx, int nx, T miny, T maxy, int ny, string prefix = ""):
        Printable(" " + prefix),
        mMinX(minx), mMaxX(maxx), mNX(nx), mMinY(miny), mMaxY(maxy), mNY(ny)
    {
        mTitle = "2D Grid";

        create();
        if (nx == 1){
            mdX = 0;
            
        }else{
            mdX = mX(1) - mX(0);
        }
        if (ny == 1){
            mdY = 0;
        }else{
            mdY = mY(1) - mY(0);
        }
    }

    Grid2D(T min, T max, T d, string prefix = ""):Printable(" " + prefix),
        mMinX(min), mMaxX(max), mdX(d), mMinY(min), mMaxY(max), mdY(d)
    {
        mTitle = "2D Grid";

        if (d == 0){
            mNX = 1;
            mNY = 1;
        }else{
            mNX = abs((mMaxX - mMinX)/mdX)+1;
            mNY = abs((mMaxY - mMinY)/mdY)+1;
        }

        create();
        mMaxX = mX(mNX-1);
        mMaxY = mY(mNY-1);
    }

    Grid2D(T min = 0, T max = 0, int n = 1, string prefix = ""):
        Printable(" " + prefix),
        mMinX(min), mMaxX(max), mNX(n), mMinY(min), mMaxY(max), mNY(n)
    {
        mTitle = "2D Grid";

        create();
        if (n == 1){
            mdX = 0;
            mdY = 0;
        }else{
            mdX = mX(1,0) - mX(0,0);
            mdY = mY(0,1) - mY(0,0);
        }
    }

    

    Col<T> X(){
        Mat<T> x = vectorise(mX);
        return x;
    }

    Col<T> Y(){
        Mat<T> y = vectorise(mY);
        return y;
    }

    T X(int i){
        return X(i);
    }

    T Y(int i){
        return Y(i);
    }
    
    svec2 operator()(int i, int j){
        svec2 xy;
        xy(0) = mX(i,j);
        xy(1) = mY(i,j);
        
        return xy;
    }
    
    T   minx() const { return mMinX; };
    T   maxx() const { return mMaxX; };
    T   dx()   const { return mdX; };
    int Nx()   const { return mNX; };
    T   miny() const { return mMinY; };
    T   maxy() const { return mMaxY; };
    T   dy()   const { return mdY; };
    int Ny()   const { return mNY; };
    
protected:
    virtual void create(){
        Col<T> x = linspace<Col<T> >(mMinX, mMaxX, mNX);
        Col<T> y = linspace<Col<T> >(mMinY, mMaxY, mNY);
        
        meshgrid<T>(mX, mY, x, y);         
    }
 
};

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

// typedefs
typedef Grid1D<double> BiasGrid;
typedef Grid1D<double> Egrid;
typedef Grid1D<double> VecGrid;
typedef Grid2D<double> MatGrid;


}
#endif	/* GRID_HPP */

