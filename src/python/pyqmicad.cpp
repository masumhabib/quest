/* 
 * File:   pyqmicad.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */

#include "pyqmicad.h"
#include "../python/boostpython.hpp"
#include "../include/qmicad.hpp"
#include "../utils/vout.h"




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

/**
 * Prints welcome message.
 */
char const* greet(){   
    static string msg;
    msg  = "      QMICAD: Quantum Mechanics Inspired Computer Aided Design";
    msg += "\n                             v" + qmicad::version + "\n";
    msg += " ------------------------------------------------------------------";
    return msg.c_str();
}

/**
 * Sets application verbosity level.
 */

void setVerbosity(int verb){
    stds::vout.appVerbosity(verbosity(verb));
}

BOOST_PYTHON_MODULE(qmicad)
{    
    def("greet", greet);
    def("setVerbosity", setVerbosity);
   
    export_Option();
    export_Printable();
    
    export_cxmat();
    export_mat();
    export_vec();

    export_VecGrid();
    
    
    export_svec();
    export_pvec();
    export_lvec();
    export_lcoord();
    export_Atom();
    export_PeriodicTable();
    export_AtomicStruct();    

    export_Workers();

    export_point();
    export_quadrilateral();
    
    export_Timer();
    
    
    export_HamParams();
    export_ham();
    export_cxham();
    export_GrapheneTbParams();
    export_GrapheneTbHam();
    export_GrapheneKpParams();
    export_GrapheneKpHam();        
    export_TISurfKpParams();
    export_TISurfKpHam();

    export_Potential();
    export_LinearPot();
     
    export_CohRgfaParams();    
    export_NegfEloop();

    export_KPoints();

    export_BandStructParams();
    export_BandStruct();    
        
    
}

}
}