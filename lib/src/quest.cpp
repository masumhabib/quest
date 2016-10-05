/* 
 * File:   quest.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 2:50 PM
 *
 * QUantum mechanics Enabled Simulation Toolset (QUEST) library is a 
 * collection of C++ classes and functions for simulation and design of 
 * nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
 * Green Function (NEGF) formalism. This library provides a common framework 
 * for NEGF so that it can be used with any empirical tight binding, k.p model, 
 * extended Huckel method and density functional theory codes. This library also
 * includes a generic empirical tight binding model and a k.p model which can be 
 * extended for any material with known parameters.
 * 
 * Unlike most of the quantum transport tools available out there, 
 * QUEST is a *not* a predefined simulator for a predefined problem in 
 * a predefined device structure. Rather, it is a flexible toolbox 
 * that allows users to write their own simulators best-suited for their 
 * own problems. 

 * Deep down, QUEST is a C++ library that uses BLAS and LAPACK for computation.
 * Yet, QUEST provides both C++ and Python interfaces so that one can choose 
 * between the elegance of C++ and ease of Python. 
 * 
 */

#include "quest.hpp"
#include <string>

namespace quest{

using utils::stds::string;
const string version = string("v") + QUEST_VERSION_MAJOR + "." + 
                       QUEST_VERSION_MINOR + "." + QUEST_VERSION_PATCH;

/**
 * Prints welcome message.
 */
char const* greet(){   
    static string msg;
    msg += "\n **************************************************************";
    msg += "\n ***                         QUEST                          ***";
    msg += "\n ***    Quantum Mechanics Inspired Computer Aided Design    ***";
    msg += "\n ***                         " + quest::version 
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




