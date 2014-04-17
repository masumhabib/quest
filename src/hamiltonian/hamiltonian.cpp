/* 
 * File:   hamiltonian.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 6, 2013, 5:52 PM
 * 
 * Description: Tight binding calculation logic and data.
 * 
 */

#include "hamiltonian.hpp"
#include "../python/boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace hamiltonian;
using namespace atoms;
using namespace utils::stds;

/**
 * Hamiltonian parameters.
 */
void export_HamParams(){    
    class_<HamParams, bases<Printable>, shared_ptr<HamParams> >("HamiltonianParams", 
            init<const string&>())
        .def_readwrite("dtol", &HamParams::dtol)
        .def("update", &HamParams::update)
    ;
}

/**
 * Real Hamiltonian base.
 */
void (ham::*ham_generate1)(const AtomicStruct&, 
        const AtomicStruct&, uint, uint) = &ham::generate;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ham_H, H, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ham_S, S, 1, 2)
void export_ham(){    
    class_<ham, shared_ptr<ham > >("Hamiltonian", no_init)
        .def("setSize", &ham::setSize)
        .def("setSizeForNegf", &ham::setSizeForNegf)
        .def("setSizeForBand", &ham::setSizeForBand)
        .def("H0", &ham::H0)
        .def("Hl", &ham::Hl)
        .def("H", &ham::H, ham_H())  
        .def("Sl", &ham::Sl)
        .def("S0", &ham::S0)
        .def("S", &ham::S, ham_S())  
        .def("genDiagBlock", &ham::genDiagBlock)
        .def("genLowDiagBlock", &ham::genLowDiagBlock)
        .def("genNearestNeigh", &ham::genNearestNeigh)
        .def("generate", ham_generate1)
    ;
}
    
/**
 * Complex Hamiltonian base.
 */
void (cxham::*cxham_generate1)(const AtomicStruct&, 
        const AtomicStruct&, uint, uint) = &cxham::generate;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxham_H, H, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxham_S, S, 1, 2)
void export_cxham(){    
    class_<cxham, shared_ptr<cxham > >("CxHamiltonian", no_init)
        .def("setSize", &cxham::setSize)
        .def("setSizeForNegf", &cxham::setSizeForNegf)
        .def("setSizeForBand", &cxham::setSizeForBand)
        .def("H0", &cxham::H0)
        .def("Hl", &cxham::Hl)
        .def("H", &cxham::H, cxham_H())  
        .def("Sl", &cxham::Sl)
        .def("S0", &cxham::S0)
        .def("S", &cxham::S, cxham_S())  
        .def("genDiagBlock", &cxham::genDiagBlock)
        .def("genLowDiagBlock", &cxham::genLowDiagBlock)
        .def("genNearestNeigh", &cxham::genNearestNeigh)
        .def("generate", cxham_generate1)
    ;
}

}
}


