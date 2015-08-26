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
using std::shared_ptr;
using std::make_shared;
using std::priority_queue;
using std::cout;
using std::endl;
using std::get;

using maths::armadillo::mat;
using maths::armadillo::vec;
using maths::armadillo::fill;
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
typedef priority_queue<Particle::ptr, vector<Particle::ptr>, ParticleComparator> 
ElectronQueue;

class ElectronBins {
public:
    ElectronBins(int nbins):bins(nbins, 0.0), totalNumElects(0.0) {}
    void putElectron(const Particle::ptr electron, int ibin) {
        double occu = electron->getOccupation();
        bins[ibin] += occu;
        totalNumElects += occu; 
    }
    void reset() {
        for (auto &occu : bins) {
            occu = 0;
        }
        totalNumElects = 0;
    }
    vec calcTrans() {
        vec T(bins.size(), fill::zeros);
        for (int i = 0; i < bins.size(); i += 1) {
           T(i) = bins[i]/totalNumElects;
       }
       return T;
    }
    double getTotalNumElects() const { return totalNumElects; };
    ElectronBins& operator+=(const ElectronBins& rhs) {
        for (int i = 0; i < bins.size(); i += 1) {
            bins[i] += rhs.bins[i];
        }
        totalNumElects += rhs.totalNumElects;

        return *this;
    }
private:
    vector<double> bins;
    double totalNumElects;
};

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

    int getMaxNumTrajsPerElect() const { return mMaxTrajsPerElect; };
    void setMaxNumTrajsPerElect(int ntrajs) { mMaxTrajsPerElect = ntrajs; };
    int getMaxNumStepsPerTraj() const { return mMaxStepsPerTraj; };
    void setMaxNumStepsPerTraj(int nsteps) { mMaxStepsPerTraj = nsteps; };
    int getNumPointsPerCycle() const { return mPtsPerCycle; };
    void setNumPointsPerCycle(int nPtsPerCycle) { mPtsPerCycle = nPtsPerCycle; };
    int getNumEdgeSeekSteps() const { return mNdtStep; };
    void setNumEdgeSeekSteps(int nsteps) { mNdtStep = nsteps; };
    double getInjecSpacing() const { return mdl; };
    void setInjecSpacing(double dl) { mdl = dl; };
    int getNumInjecDir() const { return mNth; };
    void setNumInjecDir(int ndir) { mNth = ndir; };
    double getInjecAngleSpread() const { return mAngleSpread; };
    void setInjecAngleSpread(double angle) { mAngleSpread = angle; };
    double getFermiVelo() const { return mvF; };
    void getFermiVelo(double vF) { mvF = vF; };
    ParticleType getParticleType() const { return particleType; };
    void setParticleType(ParticleType type) { particleType = type; };
    double getTimeStep() const { return dt; };
    void setTimeStep(double dt) { isAutoDt = false; this->dt = dt; };
    double getCollectionTol() const { return 1.0-mCollectionTol; };
    void setCollectionTol(double tol) { mCollectionTol = 1.0-tol; 
        mOccupationFailTol = tol*10; };

    void setDebugLvl(unsigned long debugLevel) { debug = (debugLevel > 0); };
    unsigned long getDebugLvl() { return debug ? 1L : 0L; };

private:
    tuple<mat, TrajectoryVect> calcTran(int injCont, bool saveTraj);
    tuple<int, ElectronBins, TrajectoryVect> calcTrajOneElect(point ri, 
        double thi, bool saveTraj);
    inline int calcSingleTraj(bool saveTraj, ElectronQueue &electsQu, 
        ElectronBins &bins, Trajectory& traj);
    inline void applyPotential(Particle::ptr electron);
    inline void refreshTimeStepSize(Particle::ptr electron);
    inline bool justCrossEdge(Particle::ptr electron, point ri, point rf, 
            int iEdge);
    inline bool getCloseToEdge(Particle::ptr electron, point ri, point rf, 
            int iEdge);
    inline bool stepNearEdge(Particle::ptr electron, point& ri, point& rf, 
            int iEdge, bool doCross);
    inline bool stepNearEdge2(Particle::ptr electron, point& ri, point& rf, 
            int iEdge, bool doCross);

private:
    int mMaxStepsPerTraj = 10000;  //!< maximum number of time steps before fail.
    int mMaxTrajsPerElect = 10000; //!< maximum number of trajectories for one electron
    int mPtsPerCycle = 100; //!< number of points per cyclotron cycle.
    int mNdtStep = 10; //!< maximum number of steps for determining reflection dt.
    double mvF = 1E6/nm; //!< Fermi velocity.
    double mdl = 5.0; //!< distance between two injection points in a contact
    int mNth = 50; //!< number of random directions for each contact.
    double dt = 10/(1E6/nm); //!< default time step. 
    double maxdt = 100/(1E6/nm); //!< maximum time step size
    bool isAutoDt = true; //!< Switch for automatic dt calculation
    double mB = 0; //!< Magnetic field.
    double mV = 0; //!< Electric potential.
    double mE = 0; //!< Energy of electron.
    double mAngleSpread = pi/5; //!< spread of incident angle

    bool debug = false; //!< Prints debug message if true.

    double mReflectionTol = 1E-4;
    double mTransmissionTol = 1E-4;
    double mCollectionTol = 0.99;
    double mOccupationFailTol = 1E-2;
    double mClosenessTol = 1E-2;

    static constexpr double ETOL = 1E-6;

    Device::ptr mDev; //!< Device structure.
    ParticleType particleType = ParticleType::DiracCyclotron; //!< particle type.
};

}}

#endif



