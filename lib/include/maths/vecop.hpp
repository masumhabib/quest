/* 
 * File:   vecop.hpp
 * Author: Mirza Elahi <mirza.monzur@gmail.com>
 *
 * Created on April 16, 2015, 3:12 AM
 */


#ifndef VECOP_HPP
#define	VECOP_HPP

#include "utils/Printable.hpp"
#include "utils/std.hpp"
#include "utils/vout.h"

namespace maths{

/*
 * Vector Operation similar to matlab
 */

template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vintersection(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result);

template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vdifference(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result);

template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vunion(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result);

/* vector element-wise intersection.*/
template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vintersection(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result) {
    std::vector<T> output;
    std::set_intersection(first.begin(), first.end(),
                        second.begin(), second.end(),
                        back_inserter(output));
    result = arma::conv_to<ARMA_VECTOR_TYPE<T>>::from(output);
}

/* take union of two vectors.*/
template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vdifference(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result) {
  std::vector<T> output;
  std::set_difference(first.begin(), first.end(),
                      second.begin(), second.end(),
                      back_inserter(output));
  result = arma::conv_to<ARMA_VECTOR_TYPE<T>>::from(output);
}

/* take union of two vectors.*/
template<typename T, template <typename> class ARMA_VECTOR_TYPE>
void vunion(const ARMA_VECTOR_TYPE<T> &first, const ARMA_VECTOR_TYPE<T> &second, ARMA_VECTOR_TYPE<T> &result) {
  std::vector<T> output;
  std::set_union(first.begin(), first.end(),
                  second.begin(), second.end(),
                  back_inserter(output));
  result = arma::conv_to<ARMA_VECTOR_TYPE<T>>::from(output);
}


}
#endif	/* VECOP_HPP */

