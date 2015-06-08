/* 
 * File:   pysimulator.h
 * Author: K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on May 31, 2015, 09:03 AM
 */

#ifndef TMFSC_PYLIB_PYSIMULATOR_H
#define	TMFSC_PYLIB_PYSIMULATOR_H

#include "boostpython.hpp"
#include "simulator.h"
#include "pydevice.h"

namespace qmicad { namespace python {

using tmfsc::Simulator;
using tmfsc::Device;
using tmfsc::point;
using maths::armadillo::mat;

class PySimulator : public Simulator {
public:
    PySimulator(PyDevice &dev) : Simulator(dev){};
    mat calcTrajPy(point ri, double thi, double B, 
            double EF, double V, bool saveTraj = true);
};

void export_Simulator();

}}

#endif // TMFSC_PYLIB_PYDEVICE_H 


