/**
 * file: simulator.h
 * author: K M Masum Habib
 */

#ifndef TMFSC_LIB_SIMULATOR_H
#define TMFSC_LIB_SIMULATOR_H

#include "device.h"
#include "particle.hpp"
#include "relativisticparticle.hpp"
#include "maths/constants.h"
#include "utils/random.h"
#include <tuple>
#include <iostream>

namespace qmicad{ namespace tmfsc{
using std::tuple;
using std::make_tuple;
using std::tie;
using std::cout;
using std::endl;
using maths::armadillo::mat;
using maths::armadillo::zeros;
using maths::constants::pi;
using utils::random::genNormalDist;

typedef vector<point> Trajectory;
typedef vector<Trajectory> TrajectoryVect;

class Simulator : public Printable {
public:
    Simulator(Device &dev);

    tuple<mat, TrajectoryVect> calcTran(double B, double E, double V, 
            int injCont = 0, bool saveTraj = false);
    vector<point> calcTraj(point ri, double thi, double B, 
            double EF, double V, bool saveTraj = true);
    Trajectory calcTraj(Particle& particle, bool saveTraj);
    void setMaxNumTimeStep(int nsteps) { mNSteps = nsteps; };

private:
    inline tuple<svec, svec, double> doStep(const svec &vi, double thi, 
        const svec &ri, double dth, double dt);
    void resetElectBins();
    void collectElectron(int iCont, double n = 1);
 

private:
    Device &mDev; //!< Device structure.
    int mNSteps;  //!< maximum number of time steps before fail.
    int mPtsPerCycle; //!< number of points per cyclotron cycle.
    int mNdtStep; //!< maximum number of steps for determining reflection dt.
    double mvF; //!< Fermi velocity.
    double mdl; //!< distance between two injection points in a contact
    double mNth; //!< number of random directions for each contact.

    vector<double> mElectBins; //!< Electron bins.
    double mnElects; //!< Total electrons injected.

    static constexpr double ETOL = 1E-6;
};

}}

#endif



