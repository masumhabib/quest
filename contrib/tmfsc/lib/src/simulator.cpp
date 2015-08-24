/**
 * file: simulator.cpp
 * author: K M Masum Habib
 * co-author : Mirza Elahi
 */

#include "simulator.h"

namespace qmicad { namespace tmfsc {

Simulator::Simulator(Device::ptr dev)
:mDev(dev) {
}

tuple<mat, TrajectoryVect> Simulator::calcTran(double E, double B, double V, 
        int injCont, bool saveTraj){
    int nc = mDev->numConts();
    if (injCont < 0 || injCont >= nc){
        throw invalid_argument("Contact number out of bounds");
    }

    mE = E;
    mB = B;
    mV = V;

    return calcTran(injCont, saveTraj);
}

TrajectoryVect Simulator::calcTraj(point ri, double thi, double E, 
            double B, double V, bool saveTraj) {
    if (fabs(E-V) < ETOL) {
        throw invalid_argument("EF and V are too close.");
    }
    
    mE = E;
    mB = B;
    mV = V;

    return calcTraj(ri, thi, saveTraj);
}

tuple<mat, TrajectoryVect> Simulator::calcTran(double E, double B, 
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

    mE = E;
    mB = B;

    for(int ig = 0; ig < VG.size(); ig += 1) {
        mDev->setGatePotential(ig, VG[ig]);
    }

    return calcTran(injCont, saveTraj);
}

TrajectoryVect Simulator::calcTraj(point ri, double thi, double E, 
            double B, const vector<double>& VG, bool saveTraj) {
 
    if (VG.size() != mDev->getNumGates()) {
        throw invalid_argument(" Number of gate voltages does not match"
                "number of gates");
    }

    mE = E;
    mB = B;

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
        genNormalDist(th, mAngleSpread, 0);

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
    shared_ptr<Particle> electron;
    if (particleType == ParticleType::DiracCyclotron) {
        electron = make_shared<DiracCyclotron>(ri, vi, mE, V, mB);
    } else {
        electron = make_shared<DiracElectron>(ri, vi, mE, V, mB);
    }
    refreshTimeStepSize(electron);
    mElectsQu.push(electron);

    // loop over all the electron paths created when electrons cross 
    // a transmitting boundary
    TrajectoryVect trajs;
    int itrajs = 0;
    while(!mElectsQu.empty() && itrajs < mMaxTrajsPerElect) {
        trajs.push_back(calcTraj(saveTraj));
        itrajs += 1;
    }
    if (itrajs >= mMaxTrajsPerElect) {
        if (debug) {
            cout << "-W- Maximum number of trajectories reached" << endl;
        }
    }
 
    return trajs; 
}

Trajectory Simulator::calcTraj(bool saveTraj) {
    Particle::ptr electron = mElectsQu.front();
    point ri = electron->getPos();
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
        } else if (mDev->isAbsorbEdge(iEdge)) {
            collectElectron(*electron, mDev->edgeToContIndx(iEdge));
            break;
        } else { 
            // about to crossed an edge, find the intersection and probe 
            // how close we can get to the intersection point
            if (!getCloseToEdge(*electron, ri, rf, iEdge)){
                if (debug) {
                    cout << "-W- Could not get close to edge " << iEdge << endl;
                }
                break;
            }
            // we are now close enough, now step to the point
            electron->doStep();
            // update our location
            svec r = electron->getPos();
            if (mDev->isReflectEdge(iEdge)) {
                // we were about to cross a reflecting edge, NO WAY,
                // lets reflect back
                electron->reflect(mDev->edgeNormVect(iEdge));
                rf = r;
            } else if (mDev->isTransmitEdge(iEdge)) {
                // about to cross a transmitting edge/gate boundary, let's first
                // calculate what would be transmission probability. To do so,
                // first just cross the boundary and get the potential.
                Particle::ptr transElect = electron->clone();
                if(!justCrossEdge(*transElect, r, rf, iEdge)) {
                    if (debug) {
                        cout << "-W- Could not cross the edge " << iEdge << endl;
                    }
                    break;
                }
                transElect->doStep();
                double V1 = transElect->getPot(); // potential before edge
                applyPotential(transElect);
                double V2 = transElect->getPot(); // potential after edge

                double thf, transProb, refProb;
                double thti;
                tie(thf, thti, transProb, refProb) = mDev->calcProbab(V1, V2,
                        transElect->getVel(), transElect->getEnergy(), iEdge);
                double occu = electron->getOccupation();
                //if (debug){
                //    std::cout << "-D- V1 = " << V1 << " V2 = " << V2 
                //        << " T(E) = " << transProb 
                //        << " R(E) = " << refProb << std::endl;
                //}

                // reflect?
                if (refProb > mReflectionTol && occu > mOccupationTol) {
                    Particle::ptr refElect = electron->clone();
                    refElect->reflect(mDev->edgeNormVect(iEdge));
                    refElect->setOccupation(refProb*occu);
                    mElectsQu.push(refElect);
                }
                // transmit?
                if (transProb > mTransmissionTol && occu > mOccupationTol) {
                    transElect->setOccupation(transProb*occu);
                    transElect->rotateVel(-thti - thf);
                    refreshTimeStepSize(transElect);
                    mElectsQu.push(transElect);
                }
                rf = r;
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

inline void Simulator::refreshTimeStepSize(Particle::ptr electron){
    double mydt = dt;
    if (isAutoDt) {
        // TODO: mB, mV and mE should not be a part of this class
        // TODO: get these values directly from electron. 
        double V = electron->getPot();
        double wc = mvF*mvF*nm2*mB/(mE-V); //cyclotron frequency
        mydt = abs(2*pi/wc/mPtsPerCycle); // time step in cyclotron cycle
        if (mydt > maxdt){
            mydt = maxdt;
            if (debug) {
                std::cout << "-W-  Too big time step, using default" 
                    << std::endl;
            }
        }
    }
    electron->setTimeStep(mydt);
}


inline bool Simulator::justCrossEdge(Particle& electron, point ri, point rf, 
        int iEdge) {
    return stepNearEdge(electron, ri, rf, iEdge, true); 
}

inline bool Simulator::getCloseToEdge(Particle& electron, point ri, point rf, 
        int iEdge) {
    return stepNearEdge(electron, ri, rf, iEdge, false); 
}

inline bool Simulator::stepNearEdge(Particle& electron, point& ri, point& rf, 
        int iEdge, bool doCross) {
    int itr = 0;
    double dl = doCross ? mClosenessTol/10 : -mClosenessTol/10;
    while (itr < mNdtStep) {
        point intp = mDev->intersection(iEdge, ri, rf);
        point r = electron.stepCloseToPoint(intp, dl);
        svec dr = r - intp;

        double d = distance(r, intp);
        if (mDev->intersects(iEdge, r, rf)){
            // did not cross
            if (!doCross && d < mClosenessTol){
                return true;
            }
            electron.doStep();
            ri = r;
            //rf = electron.stepCloseToPoint(intp - dr);

        } else {
            //did cross
            if (doCross && d < mClosenessTol){
                return true;
            }
            rf = r;
            //ri = electron.stepCloseToPoint(intp - dr);
        }
        // FIXME: there might be a better way than this ...
        if (!mDev->intersects(iEdge, ri, rf)) {
            return false;
        }
        itr += 1;
    }
    return false;
}

inline bool Simulator::stepNearEdge2(Particle& electron, point& ri, point& rf, 
        int iEdge, bool doCross) {
    //point rf = intp;
    int itr = 0;
    double dl = doCross ? mClosenessTol : -mClosenessTol;
    point intp = mDev->intersection(iEdge, ri, rf);
    iEdge = doCross ? iEdge : -1;
    while (itr < mNdtStep) {
        intp = electron.stepCloseToPoint(intp, dl);
        if (mDev->intersects(ri, intp) == iEdge) {
            return true;
        }
        itr += 1;
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

