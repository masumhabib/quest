/* 
 * File:   pydevice.h
 * Author: K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on May 31, 2015, 09:03 AM
 */

#ifndef TMFSC_PYLIB_PYDEVICE_H
#define	TMFSC_PYLIB_PYDEVICE_H

#include "boostpython.hpp"
#include "device.h"

namespace qmicad { namespace python {

using tmfsc::Device;
using tmfsc::svec;

class PyDevice : public Device {
public:
    PyDevice() : Device () {};
    svec edgeUnitVectPy(int indx) { return edgeUnitVect(indx); };
};

void export_Device();

}}

#endif // TMFSC_PYLIB_PYDEVICE_H 


