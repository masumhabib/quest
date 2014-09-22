/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "PyCohRgfLoop.h"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{

PyCohRgfLoop::PyCohRgfLoop(const Workers &workers, uint nb, double kT, dcmplx ieta, 
        bool orthogonal, uint nTransNeigh, string newprefix): 
        CohRgfLoop(workers, nb, kT, ieta, orthogonal, nTransNeigh, newprefix)
{
}

void PyCohRgfLoop::H0(bp::object H0, int ib, int ineigh){
    mH0(ib, ineigh) = npy2mat<dcmplx>(H0);
}

void PyCohRgfLoop::S0(bp::object S0, int ib, int ineigh){
    mS0(ib, ineigh) = npy2mat<dcmplx>(S0);
}

void PyCohRgfLoop::Hl(bp::object Hl, int ib, int ineigh){
    mHl(ib, ineigh) = npy2mat<dcmplx>(Hl);
}

void PyCohRgfLoop::Sl(bp::object Sl, int ib, int ineigh){
    mSl(ib, ineigh) = npy2mat<dcmplx>(Sl);
}

void PyCohRgfLoop::H0(bp::object H0, int ib){
    mH0(ib) = npy2mat<dcmplx>(H0);
}

void PyCohRgfLoop::S0(bp::object S0, int ib){
    mS0(ib) = npy2mat<dcmplx>(S0);
}

void PyCohRgfLoop::Hl(bp::object Hl, int ib){
    mHl(ib) = npy2mat<dcmplx>(Hl);
}

void PyCohRgfLoop::Sl(bp::object Sl, int ib){
    mSl(ib) = npy2mat<dcmplx>(Sl);
}

void PyCohRgfLoop::V(bp::object V, int ib){
    mV(ib) = npy2col<double>(V);
}

void PyCohRgfLoop::pv0(bp::object pv0, int ib, int ineigh){
    mpv0(ib, ineigh) = npy2col<double>(pv0);
}

void PyCohRgfLoop::pvl(bp::object pvl, int ib, int ineigh){
    mpvl(ib, ineigh) = npy2col<double>(pvl);
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_enableI, enableI, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_enableDOS, enableDOS, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_enablen, enablen, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_enablep, enablep, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyCohRgfLoop_save, save, 1, 2)
void (PyCohRgfLoop::*PyCohRgfLoop_H0_1)(bp::object, int, int) = &PyCohRgfLoop::H0;
void (PyCohRgfLoop::*PyCohRgfLoop_S0_1)(bp::object, int, int) = &PyCohRgfLoop::S0;
void (PyCohRgfLoop::*PyCohRgfLoop_Hl_1)(bp::object, int, int) = &PyCohRgfLoop::Hl;
void (PyCohRgfLoop::*PyCohRgfLoop_Sl_1)(bp::object, int, int) = &PyCohRgfLoop::Sl;
void (PyCohRgfLoop::*PyCohRgfLoop_H0_2)(bp::object, int) = &PyCohRgfLoop::H0;
void (PyCohRgfLoop::*PyCohRgfLoop_S0_2)(bp::object, int) = &PyCohRgfLoop::S0;
void (PyCohRgfLoop::*PyCohRgfLoop_Hl_2)(bp::object, int) = &PyCohRgfLoop::Hl;
void (PyCohRgfLoop::*PyCohRgfLoop_Sl_2)(bp::object, int) = &PyCohRgfLoop::Sl;
void (PyCohRgfLoop::*PyCohRgfLoop_V)(bp::object, int) = &PyCohRgfLoop::V;
void (PyCohRgfLoop::*PyCohRgfLoop_pv0_1)(bp::object, int, int) = &PyCohRgfLoop::pv0;
void (PyCohRgfLoop::*PyCohRgfLoop_pvl_1)(bp::object, int, int) = &PyCohRgfLoop::pvl;
void export_CohRgfLoop(){
    // ~~~~~~~~ To avoid nasty numpy segfault ~~~~~~~
    import_array(); 
    class_<CohRgfLoop, bases<Printable>, shared_ptr<CohRgfLoop>, noncopyable >("_CohRgfLoop", 
            no_init)
    ;
    
    class_<PyCohRgfLoop, bases<CohRgfLoop>, shared_ptr<PyCohRgfLoop> >("CohRgfLoop", 
            init<const Workers&, 
            optional<uint, double, dcmplx, bool, uint, string> >())
        .def("E", &PyCohRgfLoop::E)
        .def("k", &PyCohRgfLoop::k)
        .def("mu", &PyCohRgfLoop::mu)
        .def("H0", PyCohRgfLoop_H0_1)
        .def("S0", PyCohRgfLoop_S0_1)
        .def("Hl", PyCohRgfLoop_Hl_1)
        .def("Sl", PyCohRgfLoop_Sl_1)
        .def("H0", PyCohRgfLoop_H0_2)
        .def("S0", PyCohRgfLoop_S0_2)
        .def("Hl", PyCohRgfLoop_Hl_2)
        .def("Sl", PyCohRgfLoop_Sl_2)
        .def("V", PyCohRgfLoop_V)
        .def("pv0", PyCohRgfLoop_pv0_1)
        .def("pvl", PyCohRgfLoop_pvl_1)
        .def("run", &PyCohRgfLoop::run)
        .def("save", &PyCohRgfLoop::save, PyCohRgfLoop_save())
        .def("enableTE", &PyCohRgfLoop::enableTE, PyCohRgfLoop_enableTE())
        .def("enableI", &PyCohRgfLoop::enableI, PyCohRgfLoop_enableI())
        .def("enableDOS", &PyCohRgfLoop::enableDOS, PyCohRgfLoop_enableDOS())
        .def("enablen", &PyCohRgfLoop::enablen, PyCohRgfLoop_enablen())
        .def("enablep", &PyCohRgfLoop::enablep, PyCohRgfLoop_enablep())
    ;
}

}
}



