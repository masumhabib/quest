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
    double (cxhamparams::*cxhamparams_getdtol)() const = &cxhamparams::dtol;
    void (cxhamparams::*cxhamparams_setdtol)(double) = &cxhamparams::dtol;
    //double (cxhamparams::*cxhamparams_getBz)() const = &cxhamparams::Bz;
    //void (cxhamparams::*cxhamparams_setBz)(double) = &cxhamparams::Bz;    
    //void (cxhamparams::*cxhamparams_setBz2)(double, int) = &cxhamparams::Bz;
    //void (cxhamparams::*cxhamparams_setBz)(double) = static_cast< void(cxhamparams::*) (double)>(&cxhamparams::Bz);
    //BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxhamparams_getBz, cxhamparams::Bz, 0, 0)
    //BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxhamparams_setBz, cxhamparams::Bz, 1, 2)
    class_<cxhamparams, bases<Printable>, shared_ptr<cxhamparams> >("HamiltonianParams", 
            no_init)
        //.def("setBz",  static_cast< void(cxhamparams::*) (double)>(&cxhamparams::Bz), cxhamparams_setBz())
        //.def("setBz",  cxhamparams_setBz)
        .add_property("dtol", cxhamparams_getdtol, cxhamparams_setdtol)
        //.add_property("Bz", cxhamparams_getBz, cxhamparams_setBz)
        //.add_property("Bz", static_cast< double(cxhamparams::*) ()>(&cxhamparams::Bz), cxhamparams_getBz(), static_cast< void(cxhamparams::*) (double)>(&cxhamparams::Bz), cxhamparams_setBz())
    ;
    
    def("generateHamOvl", generateHamOvl, " Generates Hamiltonian and Overlap matrices.");
    
}

}
}


