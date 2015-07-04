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
#include "potential/potential.h"
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

    Simulator(Device::ptr dev);
    tuple<mat, TrajectoryVect> calcTran(double B, double E, double V, 
            int injCont = 0, bool saveTraj = false);
    vector<point> calcTraj(point ri, double thi, double B, 
            double E, double V, bool saveTraj = true);
    tuple<mat, TrajectoryVect> calcTran(double B, double E, 
            const vector<double>& VG, int injCont = 0, bool saveTraj = false);
    vector<point> calcTraj(point ri, double thi, double B, 
            double E, const vector<double>& VG, bool saveTraj = true);

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
    double getTimeStep() const { return dt; };
    void setTimeStep(double dt) { isAutoDt = false; this->dt = dt; };

private:
    tuple<mat, TrajectoryVect> calcTran(int injCont, bool saveTraj);
    Trajectory calcTraj(point ri, double thi, bool saveTraj);
    Trajectory calcTraj(Particle& particle, bool saveTraj);
    inline bool getCloseToEdge(Particle& particle, point& rf, const point& ri, 
            const point& intp);
    void resetElectBins();
    void collectElectron(int iCont, double n = 1);

private:
    Device::ptr mDev; //!< Device structure.
    int mMaxStepsPerTraj = 10000;  //!< maximum number of time steps before fail.
    int mPtsPerCycle = 100; //!< number of points per cyclotron cycle.
    int mNdtStep = 1000; //!< maximum number of steps for determining reflection dt.
    double mvF = 1E6/nm; //!< Fermi velocity.
    double mdl = 5.0; //!< distance between two injection points in a contact
    int mNth = 50; //!< number of random directions for each contact.
    double dt = 1/(1E6/nm); //!< default time step.
    bool isAutoDt = false; //!< Switch for automatic dt calculation
    double mB = 0; //!< Magnetic field.
    double mV = 0; //!< Electric potential.
    double mE = 0; //!< Energy of electron.

    vector<double> mElectBins; //!< Electron bins.
    double mnElects; //!< Total electrons injected.

    ParticleType particleType = ParticleType::DiracCyclotron; //!< particle type.
    static constexpr double ETOL = 1E-6;
};

}}

#endif



