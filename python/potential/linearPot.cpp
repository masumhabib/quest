/* 
 * File:   linearPot.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 27, 2014, 5:19 PM
 */

#include "potential/linearPot.h"
#include "python/boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace potential;

/**
 * Linear potential
 */  
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LinearPot_VLR, VLR, 3, 5)
void export_LinearPot(){
    class_<LinearPot, bases<Potential>, shared_ptr<LinearPot> >("LinearPot", 
            init<optional<AtomicStruct::ptr, const string&> >())
        .enable_pickling()
        .def("addLinearRegion", &LinearPot::addLinearRegion)
        .add_property("NLR", &LinearPot::NLR) 
        .def("VLR", &LinearPot::VLR, LinearPot_VLR()) 
    ;
}

}
}
    
    

