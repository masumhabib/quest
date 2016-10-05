/* 
 * File:   boostpython.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on March 9, 2014, 3:16 PM
 */

#ifndef BOOSTPYTHON_HPP
#define	BOOSTPYTHON_HPP

#include <Python.h>
#include <boost/python.hpp>
#include <boost/format.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/detail/wrap_python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/list.hpp>
#include <boost/mpi/python.hpp>
#include <memory>

namespace quest{
namespace python{

using namespace boost::python;
using std::shared_ptr;
using boost::noncopyable;

}
}

#endif	/* BOOSTPYTHON_HPP */

