/* 
 * File:   pyqmicad.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */

#include "pyqmicad.h"
#include "../python/boostpython.hpp"
#include "../qmicad/qmicad.hpp"



namespace qmicad{
namespace python{

using boost::mpi::communicator;
using namespace utils;
using namespace atoms;
using namespace hamiltonian;
using namespace kpoints;
using namespace potential;
using namespace band;
using namespace negf;
using namespace maths::geometry;


BOOST_PYTHON_MODULE(_utils)
{
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

BOOST_PYTHON_MODULE(_atoms)
{     
    export_svec();
    export_pvec();
    export_lvec();
    export_lcoord();
    export_Atom();
    export_PeriodicTable();
    export_AtomicStruct();     
}

BOOST_PYTHON_MODULE(_kpoints)
{
    export_KPoints();        
}

BOOST_PYTHON_MODULE(_potential)
{
    export_Potential();
    export_LinearPot();
}

BOOST_PYTHON_MODULE(_hamiltonian)
{
    export_HamParams();
    export_ham();
    export_cxham();
    export_GrapheneTbParams();
    export_GrapheneTbHam();
    export_GrapheneKpParams();
    export_GrapheneKpHam();        
    export_TISurfKpParams();
    export_TISurfKpHam();    
    export_TISurfKpParams4();
    export_TISurfKpHam4();    

}

BOOST_PYTHON_MODULE(_band)
{
    export_BandStructParams();
    export_BandStruct();    
}

BOOST_PYTHON_MODULE(_negf)
{
    export_CohRgfaParams();    
    export_NegfEloop();    
}

BOOST_PYTHON_MODULE(_qmicad)
{ 
    scope().attr("version") = version;
    def("greet", greet, " Shows the QMICAD banner.");
    def("setVerbosity", setVerbosity, " Sets the verbosity level of C++ code.");      

}

}
}