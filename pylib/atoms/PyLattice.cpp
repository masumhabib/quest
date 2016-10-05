/* 
 * File:   Lattice.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 7, 2013, 9:46 AM
 * 
 * Description: Lattice vector and lattice coordinates.
 * 
 */

#include "boostpython.hpp"
#include "atoms/Lattice.h"


/**
 * Python exporter.
 */
namespace quest{
namespace python{

void export_lvec(){
    using namespace atoms;
    
    /**
     * Lattice vector.
     */
    class_<lvec, bases<Printable>, shared_ptr<lvec> >("LVec", 
            init<optional<const string&> >())
        .def_readwrite("a1", &lvec::a1)
        .def_readwrite("a2", &lvec::a2)
        .def_readwrite("a3", &lvec::a3)
        .def("la1", &lvec::la1)
        .def("la2", &lvec::la2)
        .def("la3", &lvec::la3)
        .def(self * lcoord())
    ;
}


void export_lcoord(){
    using namespace atoms;
    /**
     * Lattice coordinate.
     */
    class_<lcoord, bases<Printable>, shared_ptr<lcoord> >("LCoord", 
            init<int, int, int, optional<const string&> >())
        .def_readwrite("n1", &lcoord::n1)
        .def_readwrite("n2", &lcoord::n2)
        .def_readwrite("n3", &lcoord::n3)
    ;
}


}
}


