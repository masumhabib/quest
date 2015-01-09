/* 
 * File:   Workers.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 16, 2014, 6:04 PM
 */

#include "parallel/Workers.h"
#include "boostpython.hpp"


namespace qmicad{namespace python{
/**
 * MPI communicator wrapper.
 */
void export_Workers(){
    using namespace parallel;
    
   class_<Workers, shared_ptr<Workers>, noncopyable>("Workers", 
          init<>())
        .def(init<int, char**>())
        .def(init<const vector<string>&>())
        .def(init<const communicator&>())
        .def("MyId", &Workers::MyId)
        .def("MasterId", &Workers::MasterId)
        .def("N", &Workers::N)
        .def("AmIMaster", &Workers::AmIMaster)
        .def("IAmMaster", &Workers::IAmMaster)
    ;     
}

}}

