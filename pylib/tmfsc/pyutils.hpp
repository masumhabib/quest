/* 
 * File:   pytutils.hpp
 * Copyright (C) 2015  K M Masum Habib <masum.habib@gmail.com>
 * 
 */

#ifndef TMFSC_PYLIB_PYUTILS_H
#define	TMFSC_PYLIB_PYUTILS_H

#include "boostpython.hpp"
#include <vector>

namespace quest{
namespace python{

using std::vector;

template <class T>
vector<T> list2vect(const list& l);

template <class T>
vector<T> list2vect(const list& l) {
    return vector<T> (stl_input_iterator<T>(l), stl_input_iterator<T>());
}


}}

#endif

