/* 
 * File:   PyBandStruct.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 16, 2014, 7:41 PM
 */

#ifndef PYBANDSTRUCT_H
#define	PYBANDSTRUCT_H

#include "../band/BandStruct.h"
#include "PyWorkers.h"

#include <python2.6/Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>


namespace qmicad{
namespace python{
using namespace boost::python;
namespace armadillo = maths::armadillo;
using namespace utils::stds;


struct PyBandStructParams: public BandStructParams{
    
};


class PyBandStruct: public BandStruct{
public:
    PyBandStruct(shared_ptr<mat> pk, const PyBandStructParams &bp, 
            const PyWorkers &workers):BandStruct(pk, bp, workers){
        
    }

};

}
}
#endif	/* PYBANDSTRUCT_H */

