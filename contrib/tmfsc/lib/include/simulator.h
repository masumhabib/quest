/**
 * file: simulator.h
 * author: K M Masum Habib
 */

#ifndef TMFSC_LIB_SIMULATOR_H
#define TMFSC_LIB_SIMULATOR_H

#include "device.h"
#include "particle.hpp"
#include "DiracElectron.hpp"
#include "DiracCyclotron.hpp"
#include "maths/constants.h"
#include "utils/random.h"
#include <tuple>
#include <iostream>
#include <memory>

namespace qmicad{ namespace tmfsc{
using std::tuple;
using std::make_tuple;
using std::tie;
using std::make_shared;
using std::shared_ptr;
using maths::armadillo::mat;
using maths::armadillo::zeros;
using maths::constants::pi;
using utils::random::genNormalDist;

typedef vector<point> Trajectory;
typedef vector<Trajectory> TrajectoryVect;

class Simulator : public Printable {
public:
    enum class ParticleType {DiracCyclotron = 0, DiracElectron = 1 };

    Simulator(Device &dev);
    tuple<mat, TrajectoryVect> calcTran(double B, double E, double V, 
            int injCont = 0, bool saveTraj = false);
    vector<point> calcTraj(point ri, double thi, double B, 
            double EF, double V, bool saveTraj = true);

    int getMaxNumStepsPerTraj() const { return mMaxStepsPerTraj; };
    void setMaxNumStepsPerTraj(int nsteps) { mMaxStepsPerTraj = nsteps; };
    int getNumPointsPerCycle() const { return mPtsPerCycle; };
    void setNumPointsPerCycle(int nPtsPerCycle) { mPtsPerCycle = nPtsPerCycle; };
    int getNumTimeStepsEdge() const { return mNdtStep; };
    void setNumTimeStepsEdge(int nsteps) { mNdtStep = nsteps; };
    double getInjecSpacing() const { return mdl; };
    void setInjecSpacing(double dl) { mdl = dl; };
    int getNumInjecDir() const { return mNth; };
    void setNumInjecDir(int ndir) { mNth = ndir; };
    double getFermiVelo() const { return mvF; };
    void getFermiVelo(double vF) { mvF = vF; };
    ParticleType getParticleType() const { return particleType; };
    void setParticleType(ParticleType type) { particleType = type; };

private:
    Trajectory calcTraj(Particle& particle, bool saveTraj);
    void resetElectBins();
    void collectElectron(int iCont, double n = 1);

private:
    Device &mDev; //!< Device structure.
    int mMaxStepsPerTraj;  //!< maximum number of time steps before fail.
    int mPtsPerCycle; //!< number of points per cyclotron cycle.
    int mNdtStep; //!< maximum number of steps for determining reflection dt.
    double mvF; //!< Fermi velocity.
    double mdl; //!< distance between two injection points in a contact
    int mNth; //!< number of random directions for each contact.

    vector<double> mElectBins; //!< Electron bins.
    double mnElects; //!< Total electrons injected.

    ParticleType particleType; //!< particle type.

    static constexpr double ETOL = 1E-6;
};

}}

#endif



