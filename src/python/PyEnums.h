/* 
 * File:   PyEnums.h
 * Author: masum
 *
 * Created on February 15, 2014, 10:34 AM
 */

#ifndef PYENUMS_H
#define	PYENUMS_H

#include "../utils/myenums.hpp"

#include <python2.6/Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>

namespace qmicad{
namespace python{
using namespace myenums;

using boost::python::pickle_suite;
using boost::python::tuple;
using namespace boost::python;

//struct OptionPickler: public pickle_suite{
//    static tuple getinitargs(const Option& o){
//        return make_tuple(o);
//    }
//}; 

}
}


#endif	/* PYENUMS_H */

