/**
 * @author K M Masum Habib <masum.habib@gmail.com>
 */

#ifndef TMFSC_TMFSC_H
#define TMFSC_TMFSC_H

#include "maths/arma.hpp"
#include "maths/svec.h"

namespace qmicad { namespace tmfsc {
typedef maths::spvec::svec svec;
typedef maths::spvec::svec point;
static constexpr double nm = 1E-9;   // nanometer
static constexpr double nm2 = nm*nm;   // nanometer squared
static constexpr double AA = 1E-10;  //  Angstrom



}}

#endif



