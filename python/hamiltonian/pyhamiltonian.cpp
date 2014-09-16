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

// Helper functions just to make boost::python happy.
double cxhamparams_getBz2(const cxhamparams &self){
    return self.Bz();
}
void cxhamparams_setBz2(cxhamparams &self, double Bz){
    self.Bz(Bz);
}

/**
 * Hamiltonian parameters.
 */
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxhamparams_setBz, Bz, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxhamparams_getBz, Bz, 0, 0)
double (cxhamparams::*cxhamparams_getdtol)() const = &cxhamparams::dtol;
void (cxhamparams::*cxhamparams_setdtol)(double) = &cxhamparams::dtol;
const PeriodicTable& (cxhamparams::*cxhamparams_getptable)() const = &cxhamparams::periodicTable;
void (cxhamparams::*cxhamparams_setptable)(const PeriodicTable&) = &cxhamparams::periodicTable;
void export_cxhamparams(){    

    
    class_<cxhamparams, bases<Printable>, shared_ptr<cxhamparams> >("HamiltonianParams", 
            no_init)
        .add_property("dtol", cxhamparams_getdtol, cxhamparams_setdtol)
        .add_property("ptable", make_function(cxhamparams_getptable, return_value_policy<reference_existing_object>()), cxhamparams_setptable)
        .add_property("Bz", &cxhamparams_getBz2, &cxhamparams_setBz2)
        .def("getBz",  static_cast< double(cxhamparams::*) () const >(&cxhamparams::Bz), cxhamparams_getBz())
        .def("setBz",  static_cast< void(cxhamparams::*) (double, int)>(&cxhamparams::Bz), cxhamparams_setBz())
    ;
    
    def("generateHamOvl", generateHamOvl, " Generates Hamiltonian and Overlap matrices.");    
}

}
}


