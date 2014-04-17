/* 
 * File:   grid.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 26, 2014, 10:08 AM
 */

#include "grid.hpp"
#include "../python/boostpython.hpp"


/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace utils;

/**
 * Vector grid
 */
void export_VecGrid(){
    class_<VecGrid, bases<Printable>, shared_ptr<VecGrid> >("VecGrid", 
            init<double, double, double, optional<const string&> >())
        .def(init<optional<double, double, int, const string&> >())
        .def("V", &VecGrid::V)        
        .def("min", &VecGrid::min)
        .def("max", &VecGrid::max)
        .def("del", &VecGrid::del)
        .def("N", &VecGrid::N)
    ;
}
}
}
    
    
    
