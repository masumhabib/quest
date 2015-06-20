/**
 * file: simulator.h
 * author: K M Masum Habib
 */

#ifndef TMFSC_LIB_SIMULATOR_H
#define TMFSC_LIB_SIMULATOR_H

#include "device.h"
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

class Simulator : public Printable {
public:
    Simulator(Device &dev);

    mat calcTran(double B, double E, double V, int injCont = 0);
    vector<point> calcTraj(point ri, double thi, double B, 
            double EF, double V, bool saveTraj = true);
    void setMaxNumTimeStep(int nsteps) { mNSteps = nsteps; };
private:
    inline tuple<svec, svec, double> doStep(const svec &vi, double thi, 
        const svec &ri, double dth, double dt);
    void resetElectCounts();
    void putElectron(int indx, int n = 0);
 

private:
    Device &mDev; //!< Device structure.
    int mNSteps;  //!< maximum number of time steps before fail.
    int mPtsPerCycle; //!< number of points per cyclotron cycle.
    int mNdtStep; //!< maximum number of steps for determining reflection dt.
    double mvF; //!< Fermi velocity.
    double mdl; //!< distance between two injection points in a contact
    double mNth; //!< number of random directions for each contact.
    bool mSaveTraj; //!< save trajectory?
    bool mShowProgress; //!< show calculation progress?

    vector<int> mElects; //!< Collected electrons.
    vector<vector<point> > mTraj; //!< Saved trajectory.

    static constexpr double ETOL = 1E-6;
};

}}

#endif



