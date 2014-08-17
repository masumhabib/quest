/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "negf/CohRgfLoop.h"
#include "negf/RgfResult.h"
#include "python/boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace negf;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableI, enableI, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableDOS, enableDOS, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enablen, enablen, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enablep, enablep, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_save, save, 1, 2)
void (CohRgfLoop::*CohRgfLoop_H0_1)(bp::object, int, int) = &CohRgfLoop::H0;
void (CohRgfLoop::*CohRgfLoop_S0_1)(bp::object, int, int) = &CohRgfLoop::S0;
void (CohRgfLoop::*CohRgfLoop_Hl_1)(bp::object, int, int) = &CohRgfLoop::Hl;
void (CohRgfLoop::*CohRgfLoop_Sl_1)(bp::object, int, int) = &CohRgfLoop::Sl;
void (CohRgfLoop::*CohRgfLoop_H0_2)(bp::object, int) = &CohRgfLoop::H0;
void (CohRgfLoop::*CohRgfLoop_S0_2)(bp::object, int) = &CohRgfLoop::S0;
void (CohRgfLoop::*CohRgfLoop_Hl_2)(bp::object, int) = &CohRgfLoop::Hl;
void (CohRgfLoop::*CohRgfLoop_Sl_2)(bp::object, int) = &CohRgfLoop::Sl;
void (CohRgfLoop::*CohRgfLoop_V)(bp::object, int) = &CohRgfLoop::V;
void export_CohRgfLoop(){
    // ~~~~~~~~ To avoid nasty numpy segfault ~~~~~~~
    import_array(); 
    
    class_<CohRgfLoop, bases<Printable>, shared_ptr<CohRgfLoop> >("CohRgfLoop", 
            init<const Workers&, 
            optional<uint, double, dcmplx, bool, string> >())
        .def("E", &CohRgfLoop::E)
        .def("k", &CohRgfLoop::k)
        .def("mu", &CohRgfLoop::mu)
        .def("H0", CohRgfLoop_H0_1)
        .def("S0", CohRgfLoop_S0_1)
        .def("Hl", CohRgfLoop_Hl_1)
        .def("Sl", CohRgfLoop_Sl_1)
        .def("H0", CohRgfLoop_H0_2)
        .def("S0", CohRgfLoop_S0_2)
        .def("Hl", CohRgfLoop_Hl_2)
        .def("Sl", CohRgfLoop_Sl_2)
        .def("V", CohRgfLoop_V)
        .def("run", &CohRgfLoop::run)
        .def("save", &CohRgfLoop::save, CohRgfLoop_save())
        .def("enableTE", &CohRgfLoop::enableTE, CohRgfLoop_enableTE())
        .def("enableI", &CohRgfLoop::enableI, CohRgfLoop_enableI())
        .def("enableDOS", &CohRgfLoop::enableDOS, CohRgfLoop_enableDOS())
        .def("enablen", &CohRgfLoop::enablen, CohRgfLoop_enablen())
        .def("enablep", &CohRgfLoop::enablep, CohRgfLoop_enablep())
    ;
}

}
}



