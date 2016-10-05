/**
 * Converts std::vector<T> to python list
 * @author Masum Habib <masum.habib@gmail.com>
 */

#ifndef PYLIB_UTILS_VECTOR2LIST_H
#define PYLIB_UTILS_VECTOR2LIST_H

#include "boostpython.hpp"
#include <vector>

namespace quest { namespace python {

//shamelessly copied from: http://stackoverflow.com/questions/6157409/stdvector-to-boostpythonlist
using std::vector;
template<class T>
list vector2list(const std::vector<T>& v){
    list l;
    typename vector<T>::const_iterator it;
    for (it = v.begin(); it != v.end(); ++it){
        l.append(*it);   
        }
    return l;  
}


}}

#endif



