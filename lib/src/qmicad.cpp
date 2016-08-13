/* 
 * File:   qmicad.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 2:50 PM
 *
 * Quantum Mechanics Inspired Computer Aided Design (QMICAD) library is a 
 * collection of C++ classes and functions for simulation and design of 
 * nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
 * Green Function (NEGF) formalism. This library provides a common framework 
 * for NEGF so that it can be used with any empirical tight binding, k.p model, 
 * extended Huckel method and density functional theory codes. This library also
 * includes a generic empirical tight binding model and a k.p model which can be 
 * extended for any material with known parameters.
 * 
 * Unlike most of the quantum transport tools available out there, 
 * QMICAD is a *not* a predefined simulator for a predefined problem in 
 * a predefined device structure. Rather, it is a flexible toolbox 
 * that allows users to write their own simulators best-suited for their 
 * own problems. 

 * Deep down, QMICAD is a C++ library that uses BLAS and LAPACK for computation.
 * Yet, QMICAD provides both C++ and Python interfaces so that one can choose 
 * between the elegance of C++ and ease of Python. 
 * 
 */

#include "qmicad.hpp"
#include <string>

namespace qmicad{

using utils::stds::string;
const string version = string("v") + QMICAD_VERSION_MAJOR + "." + 
                       QMICAD_VERSION_MINOR + "." + QMICAD_VERSION_PATCH;

/**
 * Prints welcome message.
 */
char const* greet(){   
    static string msg;
    msg += "\n **************************************************************";
    msg += "\n ***                         QMICAD                         ***";
    msg += "\n ***    Quantum Mechanics Inspired Computer Aided Design    ***";
    msg += "\n ***                         " + qmicad::version 
                                               + "                        ***";
    msg += "\n **************************************************************";
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




