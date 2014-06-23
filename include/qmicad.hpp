/* 
 * File:   qmicad.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
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




#ifndef QMICAD_HPP
#define QMICAD_HPP

#include "utils/std.hpp"
#include "utils/vout.h"
#include "utils/myenums.hpp"
#include "utils/serialize.hpp"
#include "utils/Printable.hpp"
#include "utils/ConsoleProgressBar.h"
#include "utils/Exception.h"
#include "utils/Timer.h"
#include "utils/NullDeleter.hpp"

#include "string/stringutils.h"

#include "grid/grid.hpp"

#include "maths/svec.h"
#include "maths/constants.h"
#include "maths/trace.hpp"
#include "maths/fermi.hpp"
#include "maths/geometry.hpp"

#include "atoms/Lattice.h"
#include "atoms/AtomicStruct.h"

#include "parallel/parloop.h"
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

#include "negf/computegs.h"
#include "negf/CohRgfa.h"
#include "negf/NegfEloop.h"

#include "band/BandStruct.h"

#include "armanpy/armanpy.h"

#include "config.h"

namespace qmicad{

using utils::stds::string;
extern const string version;

char const* greet();
void setVerbosity(int verb);

}

#endif /* QMICAD_HPP */
