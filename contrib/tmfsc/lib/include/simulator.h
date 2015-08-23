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
#include <queue>

namespace qmicad{ namespace tmfsc{
using std::tuple;
using std::make_tuple;
using std::tie;
using std::make_shared;
using std::shared_ptr;
using std::queue;
using std::cout;
using std::endl;
using maths::armadillo::mat;
using maths::armadillo::zeros;
using maths::constants::pi;
using utils::random::genNormalDist;

using maths::constants::pi;
using maths::armadillo::dcmplx;

typedef vector<point> Path;
struct Trajectory {
    Path path;
    double occupation = 1.0;
};
typedef vector<Trajectory> TrajectoryVect;

class Simulator : public Printable {
public:
    enum class ParticleType {DiracCyclotron = 0, DiracElectron = 1 };

    Simulator(Device::ptr dev);
    tuple<mat, TrajectoryVect> calcTran(double E, double B, double V, 
            int injCont = 0, bool saveTraj = false);
    TrajectoryVect calcTraj(point ri, double thi, double E, 
            double B, double V, bool saveTraj = true);
    tuple<mat, TrajectoryVect> calcTran(double E, double B, 
            const vector<double>& VG, int injCont = 0, bool saveTraj = false);
    TrajectoryVect calcTraj(point ri, double thi, double E, 
            double B, const vector<double>& VG, bool saveTraj = true);

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
    double getOccupationTol() const { return mOccupationTol; };
    void setOccupationTol(double tol) { mOccupationTol = tol; };

    void setDebugLvl(unsigned long debugLevel) { debug = (debugLevel > 0); };
    unsigned long getDebugLvl() { return debug ? 1L : 0L; };

private:
    tuple<mat, TrajectoryVect> calcTran(int injCont, bool saveTraj);
    TrajectoryVect calcTraj(point ri, double thi, bool saveTraj);
    Trajectory calcTraj(bool saveTraj);
    inline void applyPotential(Particle::ptr electron);
    inline void refreshTimeStepSize(Particle::ptr electron);
    inline bool justCrossEdge(Particle& electron, point ri, point rf, 
            int iEdge);
    inline bool getCloseToEdge(Particle& electron, point ri, point rf, 
            int iEdge);
    inline bool stepNearEdge(Particle& electron, point& ri, point& rf, 
            int iEdge, bool doCross);
    inline bool stepNearEdge2(Particle& electron, point& ri, point& rf, 
            int iEdge, bool doCross);
    void resetElectBins();
    void collectElectron(const Particle &electron, int iCont);

private:
    int mMaxStepsPerTraj = 10000;  //!< maximum number of time steps before fail.
    int mPtsPerCycle = 100; //!< number of points per cyclotron cycle.
    int mNdtStep = 5; //!< maximum number of steps for determining reflection dt.
    double mvF = 1E6/nm; //!< Fermi velocity.
    double mdl = 5.0; //!< distance between two injection points in a contact
    int mNth = 50; //!< number of random directions for each contact.
    double dt = 10/(1E6/nm); //!< default time step. 
    double maxdt = 100/(1E6/nm); //!< maximum time step size
    bool isAutoDt = true; //!< Switch for automatic dt calculation
    double mB = 0; //!< Magnetic field.
    double mV = 0; //!< Electric potential.
    double mE = 0; //!< Energy of electron.

    bool debug = false; //!< Prints debug message if true.

    double mReflectionTol = 1E-4;
    double mTransmissionTol = 1E-4;
    double mOccupationTol = 1E-9;
    double mClosenessTol = 1E-2;

    static constexpr double ETOL = 1E-6;

    Device::ptr mDev; //!< Device structure.
    vector<double> mElectBins; //!< Electron bins.
    double mnElects; //!< Total electrons injected.
    queue<Particle::ptr> mElectsQu; //!< Deck of electrons.
    ParticleType particleType = ParticleType::DiracCyclotron; //!< particle type.
};

}}

#endif



