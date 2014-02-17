/* 
 * File:   PyNegfEloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 8, 2014, 12:15 AM
 * 
 * Python wrapper for NegfEloop class
 * 
 */

#ifndef PYNEGFELOOP_H
#define	PYNEGFELOOP_H

#include "PyVecGrid.h"
#include "PyNegfParams.h"
#include "PyWorkers.h"
#include "../negf/NegfEloop.h"

#include <python2.6/Python.h>
#include <boost/python.hpp>
#include <boost/python/wrapper.hpp>

namespace qmicad{
namespace python{
using namespace boost::python;
/**
 * Python wrapper for NegfEloop class.
 */
class PyNegfEloop:public NegfEloop, public wrapper<NegfEloop>{
public:
   
private:

public:
    PyNegfEloop(const PyVecGrid &E, const PyNegfParams &np, 
            const PyWorkers &workers, bool saveAscii = true):
            NegfEloop(E, np, workers, saveAscii), wrapper<NegfEloop>()
    {
    };
    
    virtual void run(){
        //if (override f = this->get_override("run")){
        //    f(); 
        //}
        NegfEloop::run();
    }    
    
};
}}
#endif	/* PYNEGFELOOP_H */

