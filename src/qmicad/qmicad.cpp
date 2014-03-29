/* 
 * File:   qmicad.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 2:50 PM
 * 
 * Quantum Mechanics Inspired Computer Aided Design (QMICAD) library is a 
 * collection of C++ classes and functions for simulation and design of 
 * nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
 * Green's Function (NEGF) formalism. This library provides a common framework 
 * for NEGF so that it can be used with any empirical tight binding, k.p model, 
 * extended Huckel method and density functional theory codes. This library also
 * includes a generic empirical tight binding model and a k.p model which can be 
 * extended for any material with known parameters.
 * 
 */

#include "qmicad.hpp"
#include <string>

namespace qmicad{

using utils::stds::string;
const string version = "0.01.2";

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
    using namespace utils::stds;
    vout.appVerbosity(verbosity(verb));
}


}




