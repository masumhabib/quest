/* 
 * File:   BandStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "band/BandStruct.h"
#include "python/boostpython.hpp"


/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace band;

void export_BandStructParams(){
    class_<BandStructParams, bases<Printable>, shared_ptr<BandStructParams>, noncopyable>("BandStructParams",
        init<uint, optional<const string&> >()) 
        .def("H", &BandStructParams::setH)
        .def("S", &BandStructParams::setS)
        .def("lc", &BandStructParams::setLattCoord)
        .def_readwrite("nb", &BandStructParams::nb)
        .def_readwrite("ne", &BandStructParams::ne)
        .def_readwrite("no", &BandStructParams::no)
        .def_readwrite("lv", &BandStructParams::lv)
        .def_readwrite("isOrthogonal", &BandStructParams::isOrthogonal)
    ;
}

void export_BandStruct(){
    class_<BandStruct, shared_ptr<BandStruct>, noncopyable>("BandStruct",
            init<shared_ptr<mat>, const BandStructParams, 
            const Workers&>())
        .def("run", &BandStruct::run)
        .def("save", &BandStruct::save)
        .def("enableEigVec", &BandStruct::enableEigVec)
    ;
}

}
}


