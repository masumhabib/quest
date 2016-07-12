/* 
 * File:   pyqmicad.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */

#include "boostpython.hpp"
#include "pyqmicad.h"
#include "qmicad.hpp"



namespace qmicad{
namespace python{

using namespace boost::python;


void export_utils(){

    //map the utils namespace to a sub-module
    // make "from qmicad.utils import <whatever>" work
    object utilsModule(handle<>(borrowed(PyImport_AddModule("qmicad.utils"))));
    // make "from qmicad import utils" work
    scope().attr("utils") = utilsModule;
    //set the current scope to the new sub-module
    scope utils_scope = utilsModule;

    export_Printable();
    
    export_Option();
    
    export_Timer();
    export_cxmat();
    export_mat();
    export_vec();  

    export_VecGrid();
    
    export_Workers();

    export_point();
    export_quadrilateral();        
    
}

void export_atoms()
{     
    object atomsModule(handle<>(borrowed(PyImport_AddModule("qmicad.atoms"))));
    scope().attr("atoms") = atomsModule;
    scope atoms_scope = atomsModule;

    export_svec();
    export_pvec();
    export_lvec();
    export_lcoord();
    export_Atom();
    export_PeriodicTable();
    export_AtomicStruct();     
}

void export_kpoints()
{
    object kpointsModule(handle<>(borrowed(PyImport_AddModule("qmicad.kpoints"))));
    scope().attr("kpoints") = kpointsModule;
    scope kpoints_scope = kpointsModule;

    export_KPoints();        
}

void export_potential()
{    
    object potentialModule(handle<>(borrowed(PyImport_AddModule("qmicad.potential"))));
    scope().attr("potential") = potentialModule;
    scope potential_scope = potentialModule;

    export_Potential();
    export_LinearPot();
}

void export_hamiltonian()
{
    object hamiltonianModule(handle<>(borrowed(PyImport_AddModule("qmicad.hamiltonian"))));
    scope().attr("hamiltonian") = hamiltonianModule;
    scope hamiltonian_scope = hamiltonianModule;

    export_cxhamparams();
    export_GrapheneTbParams();
    export_GrapheneKpParams();
    export_TISurfKpParams();
    export_TISurfKpParams4();
    export_TI3DKpParams();

}

void export_band()
{
    object bandModule(handle<>(borrowed(PyImport_AddModule("qmicad.band"))));
    scope().attr("band") = bandModule;
    scope band_scope = bandModule;

    export_BandStruct();    
}

void export_negf()
{
    object negfModule(handle<>(borrowed(PyImport_AddModule("qmicad.negf"))));
    scope().attr("negf") = negfModule;
    scope negf_scope = negfModule;

    export_CohRgfLoop();    
}

void export_tmfsc()
{ 
    using namespace qmicad::tmfsc;

    object tmfscModule(handle<>(borrowed(PyImport_AddModule("qmicad.tmfsc"))));
    scope().attr("tmfsc") = tmfscModule;
    scope tmfsc_scope = tmfscModule;

    // add global variables
    object package = scope();
    package.attr("nm") = nm;
    package.attr("AA") = AA;
    package.attr("EDGE_ABSORB") = Edge::EDGE_ABSORB;
    package.attr("EDGE_REFLECT") = Edge::EDGE_REFLECT;
    package.attr("EDGE_TRANSMIT") = Edge::EDGE_TRANSMIT;

    export_Device();
    export_Simulator();
    
}

BOOST_PYTHON_MODULE(qmicad)
{ 
    // create qmicad package
    object package = scope();
    package.attr("__path__") = "qmicad";
    package.attr("version") = version;

    def("greet", greet, " Shows the QMICAD banner.");
    def("setVerbosity", setVerbosity, " Sets the verbosity level of C++ code.");

    
    export_npyarma();    
    export_utils();
    export_atoms();
    export_kpoints();
    export_potential();
    export_hamiltonian();
    export_band();
    export_negf();
    export_tmfsc();
}

}
}
