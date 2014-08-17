/* 
 * File:   hamiltonian.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 6, 2013, 5:52 PM
 * 
 * Description: Tight binding calculation logic and data.
 * 
 */

#include "hamiltonian/hamiltonian.hpp"
#include "boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace hamiltonian;
using namespace atoms;
using namespace utils::stds;
namespace bp = boost::python;

bp::tuple generateHamOvl(const HamParams<cxmat> &p, const AtomicStruct &bi, 
        const AtomicStruct &bj)
{
    cxmat H, S;
    generateHamOvl(H, S, p, bi, bj);
    
    return bp::make_tuple(H, S);
}

/**
 * Hamiltonian parameters.
 */
void export_cxhamparams(){    
    double (cxhamparams::*cxhamparams_getdtol)() = &cxhamparams::dtol;
    void (cxhamparams::*cxhamparams_setdtol)(double) = &cxhamparams::dtol;
    class_<cxhamparams, bases<Printable>, shared_ptr<cxhamparams> >("HamiltonianParams", 
            no_init)
        .add_property("dtol", cxhamparams_getdtol, cxhamparams_setdtol)
    ;
    
    def("generateHamOvl", generateHamOvl, " Generates Hamiltonian and Overlap matrices.");
    
}

}
}


