/* 
 * File:   quest.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 2:50 PM
 * 
 * QUantum mechanics Enabled Simulation Toolset (QUEST) is a 
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
 * Yet, QUEST provides both C++ and Python interfaces so that one can choose between 
 * the elegance of C++ and ease of Python. 
 * 
 */




#ifndef QUEST_HPP
#define QUEST_HPP

#include "utils/std.hpp"
#include "utils/vout.h"
#include "utils/myenums.hpp"
#include "utils/serialize.hpp"
#include "utils/Printable.hpp"
#include "utils/ConsoleProgressBar.h"
#include "utils/Exception.h"
#include "utils/Timer.h"
#include "utils/NullDeleter.hpp"
#include "utils/stringutils.h"

#include "maths/svec.h"
#include "maths/constants.h"
#include "maths/trace.hpp"
#include "maths/fermi.hpp"
#include "maths/geometry.hpp"
#include "maths/grid.hpp"
#include "maths/linspace.hpp"

#include "atoms/Lattice.h"
#include "atoms/AtomicStruct.h"

#include "parallel/Workers.h"

#include "kpoints/KPoints.h"

#include "hamiltonian/hamiltonian.hpp"
#include "hamiltonian/tb/graphenetb.h"
#include "hamiltonian/kp/graphenekp.h"
#include "hamiltonian/kp/tikp.h"
#include "hamiltonian/kp/tikp4.h"

#include "potential/terminal.h"
#include "potential/potential.h"
#include "potential/linearPot.h"

#include "band/BandStruct.h"

#include "negf/computegs.h"
#include "negf/CohRgfa.h"
#include "negf/CohRgfLoop.h"

#include "tmfsc/device.h"
#include "tmfsc/simulator.h"

#include "config.h"

namespace quest {

using utils::stds::string;
extern const string version;

char const* greet();
void setVerbosity(int verb);

}

#endif /* QUEST_HPP */
