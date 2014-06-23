/**
 *       Filename:  armanpy.h
 *
 *    Description:  Conversion utility between armadillo and numpy.
 *
 *        Created:  06/23/2014 03:59:24 PM
 *         Author:  K M Masum Habib (KMMH), masum.habib@gmail.com
 *      Copyright:  (c) 2014 K M Masum Habib<masum.habib@gmail.com>.
 *
 */


#ifndef ARMANPY_H
#define ARMANPY_H

#include "python/boostpython.hpp"
#include "maths/arma.hpp"

#include <numpy/ndarrayobject.h>

namespace qmicad{
namespace armanpy{

    using namespace boost::python;
    using namespace maths::armadillo;

    object vec_to_numpy();


}
}
#endif



