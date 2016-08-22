/* 
 * File:   std.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 9, 2014, 12:46 AM
 */

#ifndef STD_HPP
#define	STD_HPP

#include <ios>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <map>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <cmath>
#include <memory>

namespace utils{
namespace stds {

using std::cout;
using std::cin;
using std::cerr;
using std::endl;

using std::ios;
using std::ios_base;
using std::basic_ostream;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::stringstream;
using std::ostringstream;
using std::istringstream;

using std::string;
using std::vector;
using std::list;
using std::pair;
using std::map;

using std::abs;

using std::swap;

using std::invalid_argument;
using std::runtime_error;

using std::binary_function;

using std::make_shared;
using std::shared_ptr;
using std::unique_ptr;
using std::static_pointer_cast;

template<typename NumericType>
const NumericType& max (const std::vector<NumericType>& vec);

template<typename NumericType>
const NumericType& min (const std::vector<NumericType>& vec);

template<typename NumericType>
double mean (const std::vector<NumericType>& vec);

template<typename NumericType>
double std_dev (const std::vector<NumericType>& vec);

template<typename NumericType>
std::vector<NumericType> linspace (NumericType min, NumericType max, 
        size_t count);

template<typename NumericType>
std::vector<NumericType> linspace (NumericType min, NumericType max, 
        NumericType delta);


// ======= IMPLEMENTATION ========
template<typename NumericType>
const NumericType& max (const std::vector<NumericType>& vec) {
    return *std::max_element(std::begin(vec), std::end(vec));
}

template<typename NumericType>
const NumericType& min (const std::vector<NumericType>& vec) {
    return *std::min_element(std::begin(vec), std::end(vec));
}

template<typename NumericType>
double mean (const std::vector<NumericType>& vec) {
    double sum = 0;

	for (auto n : vec) {
	    sum += n; 
    }
    
    return sum/vec.size();
}

template<typename NumericType>
double std_dev (const std::vector<NumericType>& vec) {
    double sum = 0;
    double sqrd_sum = 0;

	for (auto n : vec) {
	    sum += n; 
	    sqrd_sum += n*n;
    }
    
    return std::sqrt(sqrd_sum/vec.size() - sum*sum/vec.size()/vec.size());
}


template<typename NumericType>
std::vector<NumericType> linspace (NumericType min, NumericType max, NumericType delta) {
    if (min > max) {
        throw std::invalid_argument("linspace(): min > max.");
    }

    std::vector<NumericType> vec;
    if (delta == 0) {
        return vec;
    }

    // for performance boost
    size_t count = size_t((max - min)/delta) + 1;
    vec.reserve(count);

    NumericType x = min;
    while (x <= max) {
        vec.push_back(x);
        x += delta;
    }

    return vec;
}

template<typename NumericType>
std::vector<NumericType> linspace (NumericType min, NumericType max, 
size_t count) {
    std::vector<NumericType> vec;
    vec.reserve(count);

    if (count == 0) {
        return vec;
    //} else if (count == 1) {
    //    vec.push_back(min);
    //    return vec;
    }

    NumericType dx = (max - min)/(count-1);
    //for (int i = 0; i < count; i += 1) {
    //    NumericType x = min + i*dx; 
    //    vec.push_back(x);
    //}
    return linspace(min, max, dx);

    return vec;

}

}
}


#endif	/* STD_HPP */

