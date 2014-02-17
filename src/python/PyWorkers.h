/* 
 * File:   PyWorkers.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 16, 2014, 6:34 PM
 */

#ifndef PYWORKERS_H
#define	PYWORKERS_H

#include "../parallel/Workers.h"

#include <python2.6/Python.h>
#include <boost/python.hpp>
#include <boost/python/wrapper.hpp>

namespace qmicad{
namespace python{
using namespace boost::python;
using namespace boost::mpi;
using utils::Workers;

class PyWorkers: public Workers {
public:
    PyWorkers(const communicator &workers):Workers(workers){
    }
};


}
}
#endif	/* PYWORKERS_H */

