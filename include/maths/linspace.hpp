/* 
 * File:   linspace.hpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on October 6, 2014, 3:06 PM
 */

#ifndef LINSPACE_HPP
#define	LINSPACE_HPP

#include "utils/Printable.hpp"
#include "utils/std.hpp"
#include "utils/vout.h"
#include "maths/svec.h"
#include <limits>

namespace maths{
using maths::armadillo::Col;
using maths::armadillo::Mat;
using std::stringstream;
using std::numeric_limits;


/*
 * Generates linear grid
 */
template<class T>
Col<T> linspace(T start, T end, T delta);
template<class T>
void linspace(Col<T>& col, T start, T end, T delta);


/* Definitions */
template<class T>
Col<T> linspace(T start, T end, T delta){
    long N;
    if (delta <= numeric_limits<T>::min()){
        N = 1;
        end = start;
    }else{
        N = std::abs((long)((end - start)/delta))+1;
        end = start + (N - 1)*delta;
    }
    
    using namespace utils::stds;
    dout << "DBG: beg = " << start <<  " end = " << end << " delta = " << delta << endl;
    return maths::armadillo::linspace<Col<T> >(start, end, N);
}

template<class T>
void linspace(Col<T>& col, T start, T end, T delta){
    col = maths::armadillo::linspace<T>(start, end, delta);
}

}
#endif	/* LINSPACE_HPP */

