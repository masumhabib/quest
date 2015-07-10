/**
 * file: simulator.cpp
 * author: K M Masum Habib
 */

#include "simulator.h"

namespace qmicad { namespace tmfsc {

Simulator::Simulator(Device::ptr dev)
:mDev(dev) {
}

tuple<mat, TrajectoryVect> Simulator::calcTran(double B, double E, double V, 
        int injCont, bool saveTraj){
    int nc = mDev->numConts();
    if (injCont < 0 || injCont >= nc){
        throw invalid_argument("Contact number out of bounds");
    }

    mB = B;
    mE = E;
    mV = V;

    return calcTran(injCont, saveTraj);
}

TrajectoryVect Simulator::calcTraj(point ri, double thi, double B, 
            double E, double V, bool saveTraj) {
    if (fabs(E-V) < ETOL) {
        throw invalid_argument("EF and V are too close.");
    }
    
    mB = B;
    mE = E;
    mV = V;

    return calcTraj(ri, thi, saveTraj);
}

tuple<mat, TrajectoryVect> Simulator::calcTran(double B, double E, 
        const vector<double>& VG, int injCont, bool saveTraj)
{
    int nc = mDev->numConts();
    if (injCont < 0 || injCont >= nc){
        throw invalid_argument("Contact number out of bounds");
    }

    if (VG.size() != mDev->getNumGates()) {
        throw invalid_argument(" Number of gate voltages does not match"
                "number of gates");
    }

    mB = B;
    mE = E;

    for(int ig = 0; ig < VG.size(); ig += 1) {
        mDev->setGatePotential(ig, VG[ig]);
    }

    return calcTran(injCont, saveTraj);
}

TrajectoryVect Simulator::calcTraj(point ri, double thi, double B, 
            double E, const vector<double>& VG, bool saveTraj) {
 
    if (VG.size() != mDev->getNumGates()) {
        throw invalid_argument(" Number of gate voltages does not match"
                "number of gates");
    }

    mB = B;
    mE = E;

    for(int ig = 0; ig < VG.size(); ig += 1) {
        mDev->setGatePotential(ig, VG[ig]);
    }

    return calcTraj(ri, thi, saveTraj);
}

tuple<mat, TrajectoryVect> Simulator::calcTran(int injCont, bool saveTraj){
    vector<point> injPts = mDev->createPointsOnCont(injCont, mdl);
    double th0 = mDev->contDirctn(injCont) + pi/2;
    TrajectoryVect trajs;

    // reset the bins
    resetElectBins();
    int npts = injPts.size();

    for (int ip = 0; ip < npts; ip += 1) {
        point ri = injPts[ip];
        vector<double> th(mNth);
        genNormalDist(th, pi/5, 0);

        for (double thi:th){
            if (abs(thi) < (pi/2.0-pi/20.0)) {
                TrajectoryVect traj = calcTraj(ri, th0 + thi, saveTraj);
                if (saveTraj) {
                    trajs.insert(trajs.end(), traj.begin(), traj.end());
                }
            }
        }
    }

    int nc = mDev->numConts();
    mat TE = zeros<mat>(nc,nc);
    for (int ic = 0; ic < nc; ic += 1) {
        TE(injCont, ic) = mElectBins[ic]/mnElects;
        TE(ic, injCont) = mElectBins[ic]/mnElects;
    }

    return make_tuple(TE, trajs);
}

TrajectoryVect Simulator::calcTraj(point ri, double thi, bool saveTraj) {
    double V = mV;
    if (mDev->getNumGates() > 0) {
        V = mDev->getPotAt(ri);
    }

    svec vi = {mvF*cos(thi), mvF*sin(thi)}; // inital velocity
    if (isAutoDt) {
        double wc = mvF*nm*mvF*nm*mB/(mE-mV); //cyclotron frequency
        dt = abs(2*pi/wc/mPtsPerCycle); // time step in cyclotron cycle
        if (dt > 1E-3) {
            std::cout << "-W-  Too big time step" << std::endl;
        }
    }

    shared_ptr<Particle> electron;
    if (particleType == ParticleType::DiracCyclotron) {
        electron = make_shared<DiracCyclotron>(ri, vi, mE, V, mB);
    } else {
        electron = make_shared<DiracElectron>(ri, vi, mE, V, mB);
    }
    electron->setTimeStep(dt);
    mElectsQu.push(electron);

    // loop over all the electron paths created when electrons cross 
    // a transmitting boundary
    TrajectoryVect trajs;
    while(!mElectsQu.empty()) {
        trajs.push_back(calcTraj(saveTraj));
    }
 
    return trajs; 
}

Trajectory Simulator::calcTraj(bool saveTraj) {
    Particle::ptr electron = mElectsQu.front();
    point ri = electron->getPos();
    svec vi = electron->getVel();
    svec rf = ri;

    Trajectory traj;
    if (saveTraj) {
        traj.path.push_back(ri);
    }

    int ii = 0;
    while (ii < mMaxStepsPerTraj) {
        // get the next position we are about to take, but do not step yet
        rf = electron->nextPos();
        // check if electron crossed an edge
        int iEdge = mDev->intersects(ri, rf);
        if (iEdge == -1) { 
            // no crossing, continue
            electron->doStep();
            applyPotential(electron);
        } else { 
            // about to crossed an edge, find the intersection and probe 
            // how close we can get to the intersection point
            point intp = mDev->intersection(iEdge, ri, rf);
            if (!getCloseToEdge(*electron, ri, intp, iEdge)){
                std::cout << "-W- Could not get close to edge!" << std::endl;
                break;
            }
            // we are now close enough, now step to the point
            electron->doStep();
            // update our location
            rf = electron->getPos();
            ri = rf;

            if (mDev->isReflectEdge(iEdge)) {
                // we were about to cross a reflecting edge, NO WAY,
                // lets reflect back
                electron->reflect(mDev->edgeNormVect(iEdge));
            } else if (mDev->isTransmitEdge(iEdge)) {
                // about to cross a transmitting edge/gate boundary, let's first
                // calculate what would be transmission probability. To do so,
                // first just cross the boundary and get the potential.
                Particle::ptr transElect = electron->clone();
                if(!justCrossEdge(*transElect, ri, intp, iEdge)) {
                    std::cout << "-W- Could not cross the edge!" << std::endl;
                    break;
                }
                transElect->doStep();
                double V1 = transElect->getPot(); // potential before edge
                applyPotential(transElect);
                double V2 = transElect->getPot(); // potential after edge

                double thf, transProb, refProb;
                tie(thf, transProb, refProb) = mDev->calcProbab(V1, V2, 
                        transElect->getVel(), transElect->getEnergy(), iEdge);
                double occu = electron->getOccupation();
                // reflect?
                if (refProb > REFLECTION_TOL && occu > OCCUPATION_TOL) {
                    Particle::ptr refElect = electron;
                    refElect->reflect(mDev->edgeNormVect(iEdge));
                    refElect->setOccupation(refProb*occu);
                    mElectsQu.push(refElect);
                }
                // transmit?
                if (transProb > TRANSMISSION_TOL && occu > OCCUPATION_TOL) {
                    transElect->setOccupation(transProb*occu);
                    mElectsQu.push(transElect);
                }
                break;
            } else if (mDev->isAbsorbEdge(iEdge)) {
                collectElectron(*electron, mDev->edgeToContIndx(iEdge));
                break;
            }
        }    

        // reset ourselves, ready for the next step
        if (saveTraj) {
            traj.path.push_back(rf);
        }
        ri = rf;
        ii += 1;
    }

    // last point that we have missed.
    if (saveTraj) {
        traj.path.push_back(rf);
        traj.occupation = electron->getOccupation();
    }
 
    mElectsQu.pop();
    return traj;
}

inline void Simulator::applyPotential(Particle::ptr electron) {
    if (mDev->getNumGates() > 0 ) {
        electron->setPot(mDev->getPotAt(electron->getPos()));
    }
}

inline bool Simulator::justCrossEdge(Particle& electron, 
        const point& ri, const point& intp, int iEdge) {
    point rf = electron.stepCloseToPoint(intp, CLOSENESS_TOL);
    if (mDev->intersects(ri, rf) == iEdge) {
        return true;
    }
    return false;
}

inline bool Simulator::getCloseToEdge(Particle& electron, 
        const point& ri, const point& intp, int iEdge) 
{
    point rf = electron.stepCloseToPoint(intp, -CLOSENESS_TOL);
    if (mDev->intersects(ri, rf) != iEdge) {
        return true;
    }
    return false;
}

void Simulator::collectElectron(const Particle &electron, int iCont){
    if (iCont < mElectBins.size()) {
        mElectBins[iCont] += electron.getOccupation();
        mnElects += electron.getOccupation();
    }
}

void Simulator::resetElectBins() {
    mElectBins.resize(mDev->numConts());
    for (auto &bin : mElectBins) {
        bin = 0;
    }
    mnElects = 0;
}


}}

