/* 
 * File:   pysimulator.h
 * Author: K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on May 31, 2015, 09:03 AM
 */

#ifndef TMFSC_PYLIB_PYSIMULATOR_H
#define	TMFSC_PYLIB_PYSIMULATOR_H

#include "boostpython.hpp"
#include "pyutils.hpp"
#include "simulator.h"
#include "pydevice.h"
#include <vector>

namespace qmicad { namespace python {

using tmfsc::Simulator;
using tmfsc::Device;
using tmfsc::point;
using tmfsc::TrajectoryVect;
using tmfsc::Trajectory;
using maths::armadillo::mat;
using std::vector;

struct PyTrajectory {
    mat path;
    double occupation = 1.0;
    mat getPath() const { return path; }; // just to export mat to python as numpy array
};

class PySimulator : public Simulator {
public:
    PySimulator(Device::ptr dev) : Simulator(dev){};
    list calcTrajPy(point ri, double thi, double B, double EF, double V, 
            bool saveTraj = true);
    tuple calcTranPy(double B, double E, double V, int injCont = 0, 
            bool saveTraj = false);
    list calcTrajPy2(point ri, double thi, double B, double E, const list& VG, 
            bool saveTraj = true);
    tuple calcTranPy2(double B, double E, const list& VG, int injCont = 0, 
            bool saveTraj = false);

    int getParticleTypePy();
    void setParticleTypePy(int type);

private:
     PyTrajectory Traj2PyTraj (const Trajectory& traj);
     list TrajVect2List(const TrajectoryVect& trajs);
 
};

void export_Simulator();

}}

#endif // TMFSC_PYLIB_PYDEVICE_H 


