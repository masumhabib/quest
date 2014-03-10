/* 
 * File:   boostpython.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on March 9, 2014, 3:16 PM
 */

#ifndef BOOSTPYTHON_HPP
#define	BOOSTPYTHON_HPP

#include <python2.6/Python.h>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/detail/wrap_python.hpp>

namespace qmicad{
namespace python{

using namespace boost::python;
using boost::shared_ptr;
using boost::noncopyable;

}
}

#endif	/* BOOSTPYTHON_HPP */

