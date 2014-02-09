/* 
 * File:   qmicad.hpp
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




#ifndef QMICAD_HPP
#define QMICAD_HPP

#include <string>

namespace qmicad{
    static const std::string version = "0.1";
}

#include "../utils/myenums.hpp"
#include "../utils/serialize.hpp"
#include "../utils/Printable.hpp"
#include "../utils/ConsoleProgressBar.h"
#include "../utils/Exception.h"

#include "../string/stringutils.h"

#include "../maths/svec.h"
#include "../maths/constants.h"
#include "../maths/trace.hpp"
#include "../maths/fermi.hpp"

#include "../atoms/Lattice.h"
#include "../atoms/Atoms.h"
#include "../grid/grid.hpp"

#include "../qm/genHam.hpp"
#include "../qm/kp/genKpAtoms.h"
#include "../qm/kp/kpParams.hpp"
#include "../qm/kp/graphenekp.h"
#include "../qm/kp/tikp.h"

#include "../qm/tb/graphenetb.hpp"

#include "../potential/terminal.h"
#include "../potential/potential.h"
#include "../potential/linearPot.h"

#include "../negf/computegs.h"
#include "../negf/NEGF.h"
#include "../parallel/parloop.h"
#include "../negf/NegfEloop.h"

#include "../simulations/Device.h"

#include "../python/PyNegfEloop.h"
#include "../python/PyNegfParams.h"


#endif /* QMICAD_HPP */
