/* 
 * File:   pydevice.h
 * Author: K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on May 31, 2015, 09:03 AM
 */

#ifndef TMFSC_PYLIB_PYDEVICE_H
#define	TMFSC_PYLIB_PYDEVICE_H

#include "boostpython.hpp"
#include "tmfsc/device.h"

namespace qmicad { namespace python {

using tmfsc::Device;
using tmfsc::svec;
using tmfsc::point;
using tmfsc::Edge;

class PyDevice : public Device {
public:
    PyDevice() : Device () {};
    svec edgeUnitVectPy(int indx) { return edgeUnitVect(indx); };
    point contMidPointPy(int indx);
    point edgeMidPointPy(int indx);
};

void export_Device();

}}

#endif // TMFSC_PYLIB_PYDEVICE_H 


