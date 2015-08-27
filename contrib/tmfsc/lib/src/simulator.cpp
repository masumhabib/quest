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

    auto result = calcTrajOneElect(ri, thi, saveTraj);
    return get<2>(result);
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

    auto result = calcTrajOneElect(ri, thi, saveTraj);
    return get<2>(result);
}

tuple<mat, TrajectoryVect> Simulator::calcTran(int injCont, bool saveTraj){
    return calcTranRandom(injCont, saveTraj);
}

inline tuple<mat, TrajectoryVect> Simulator::calcTranRandom(int injCont, 
        bool saveTraj){
    int nconts = mDev->numConts();
    point r0 = mDev->contMidPoint(injCont);
    svec contVect = mDev->contUnitVect(injCont);
    double contWidth = mDev->contWidth(injCont);
    double th0 = mDev->contDirctn(injCont) + pi/2;

    TrajectoryVect trajs;
    ElectronBins electBins(nconts);
    double maxError = numeric_limits<double>::max();
    row prevTE = row(nconts);
    prevTE.fill(numeric_limits<double>::min());
    int ip = 0;

    while (maxError < mTransmissionTol) {
        double distance = getUniformRand(-contWidth/2, contWidth/2); 
        double angle = getGaussianRand(mAngleSpread, 0, -pi/2+mAngleLimit, 
                pi/2-mAngleLimit);
        svec position = r0 + distance * contVect;
        angle = th0 + angle;

        int status;
        TrajectoryVect traj;
        ElectronBins bin(nconts);
        tie(status, bin, traj) = calcTrajOneElect(position, angle, saveTraj);
        
        if (status == -1) {
            continue;
        }
        if (saveTraj) {
            trajs.insert(trajs.end(), traj.begin(), traj.end());
        }
       
        electBins += bin;
        row newTE = electBins.calcTransVec();
        maxError = max((newTE - prevTE)/prevTE);
        prevTE = newTE;

        ip += 1;
        if (ip >= mMaxNumInjPoints) {
            if (debug) {
                cout << "-W- Transmission calculation did not converge within "
                     << mTransmissionTol << " in " << ip << " injections." << endl;
            }
            break;
        }
    }

    mat TE = electBins.calcTransMat(nconts);
    return make_tuple(TE, trajs);
}

inline tuple<mat, TrajectoryVect> Simulator::calcTranSemiRandom(int injCont, 
        bool saveTraj){
    int nconts = mDev->numConts();
    vector<point> injPts = mDev->createPointsOnCont(injCont, mdl);
    int npts = injPts.size();
    double th0 = mDev->contDirctn(injCont) + pi/2;

    if (npts * mNth > mMaxNumInjPoints) {
        throw invalid_argument("Too many injection points, please reduce it.");
    }

    TrajectoryVect trajs;
    ElectronBins electBins(nconts);

    for (int ip = 0; ip < npts; ip += 1) {
        point ri = injPts[ip];
        vector<double> th(mNth);
        genNormalDist(th, mAngleSpread, 0);

        for (double thi:th){
            if (abs(thi) < (pi/2.0-mAngleLimit)) {
                int status;
                TrajectoryVect traj;
                ElectronBins bin(nconts);
                tie(status, bin, traj) = calcTrajOneElect(ri, th0 + thi, 
                        saveTraj);
                
                if (status == -1) {
                    continue;
                }
                
                electBins += bin;

                if (saveTraj) {
                    trajs.insert(trajs.end(), traj.begin(), traj.end());
                }
            }
        }
    }

    mat TE = electBins.calcTransMat(nconts);
    return make_tuple(TE, trajs);
}


inline tuple<int, ElectronBins, TrajectoryVect> Simulator::calcTrajOneElect(
        point ri, double thi, bool saveTraj) {
    int status = 0;
    double V = mV;
    if (mDev->getNumGates() > 0) {
        V = mDev->getPotAt(ri);
    }

    svec vi = {mvF*cos(thi), mvF*sin(thi)}; // inital velocity
    Particle::ptr electron;
    if (particleType == ParticleType::DiracCyclotron) {
        electron = make_shared<DiracCyclotron>(ri, vi, mE, V, mB);
    } else {
        electron = make_shared<DiracElectron>(ri, vi, mE, V, mB);
    }
    refreshTimeStepSize(electron);
    ElectronQueue electsQu;
    electsQu.push(electron);

    // loop over all the electron paths created when electrons cross 
    // a transmitting boundary and calculate the trajectory for each of these
    // electrons; loop terminates if (1) significant part of the electron 
    // has already been collected or (2) if we have reached maximum number
    // of trajectories allowed or (3) the queue is empty -- also implies (1).
    TrajectoryVect trajs;
    ElectronBins electBins(mDev->numConts());
    int itrajs = 0;
    while(!electsQu.empty()) {
        Trajectory traj;
        status = calcSingleTraj(saveTraj, electsQu, electBins, traj);
        if (status == -1) {
            break;
        }
        if (saveTraj) {
            trajs.push_back(traj);
        }
        if (electBins.getTotalNumElects() > mCollectionTol) {
            status = 0;
            break;
        }

        itrajs += 1;
        if (itrajs >= mMaxTrajsPerElect) {
            status = -2;
            if (debug) {
                cout << "-W- Maximum number of trajectories reached" << endl;
            }
            break;
        }
    }
 
    return make_tuple(status, electBins, trajs); 
}

inline int Simulator::calcSingleTraj(bool saveTraj, ElectronQueue &electsQu, 
        ElectronBins &bins, Trajectory& traj) {
    int status = 0;
    Particle::ptr electron = electsQu.top();
    electsQu.pop();

    point ri = electron->getPos();
    svec rf = ri;

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
            // crossed an absorbing edge, collect it and then we'll be done.
            bins.putElectron(electron, mDev->edgeToContIndx(iEdge));
            status = 1;
            break;
        } else { 
            // about to cross a reflecting or transmitting edge, find the 
            // intersection and probe how close we can get to the 
            // intersection point
            if (!getCloseToEdge(electron, ri, rf, iEdge)){
                if (debug) {
                    cout << "-W- Could not get close to edge " << iEdge << endl;
                }
                // we failed to get close to edge, this part of the electron
                // will be discarded; therefor, we can not continue with this
                // electron if it carries a significant portion ...
                if (electron->getOccupation() > mOccupationFailTol) {
                    status = -1;
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
                if(!justCrossEdge(transElect, r, rf, iEdge)) {
                    if (debug) {
                        cout << "-W- Could not cross the edge " << iEdge << endl;
                    }
                    // if this electron carries a significant occupation,
                    // we should no longer use it in our transmission 
                    // calculation
                    if (transElect->getOccupation() > mOccupationFailTol) {
                        status = -1;
                    }
                    break;
                }
                // cross the transmission boundary
                transElect->doStep();
                // get the potentials at the both sides of the transmitting edge
                double V1 = transElect->getPot(); // potential before edge
                applyPotential(transElect);
                double V2 = transElect->getPot(); // potential after edge

                // calculate the transmission and reflection probability
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

                // reflect? if yes, reflect it and put it in the queue
                if (refProb > mReflectionTol) {
                    Particle::ptr refElect = electron->clone();
                    refElect->reflect(mDev->edgeNormVect(iEdge));
                    refElect->setOccupation(refProb*occu);
                    electsQu.push(refElect);
                }
                // transmit? if yes, transmit it and put it in the queue
                if (transProb > mTransmissionTol) {
                    transElect->setOccupation(transProb*occu);
                    transElect->rotateVel(-thti - thf);
                    refreshTimeStepSize(transElect);
                    electsQu.push(transElect);
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

    // append the last point to trajectory that we have missed.
    if (saveTraj) {
        traj.path.push_back(rf);
        traj.occupation = electron->getOccupation();
    }
    return status;
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


inline bool Simulator::justCrossEdge(Particle::ptr electron, point ri, 
        point rf, int iEdge) {
    return stepNearEdge(electron, ri, rf, iEdge, true); 
}

inline bool Simulator::getCloseToEdge(Particle::ptr electron, point ri, 
        point rf, int iEdge) {
    return stepNearEdge(electron, ri, rf, iEdge, false); 
}

inline bool Simulator::stepNearEdge(Particle::ptr electron, point& ri, 
        point& rf, int iEdge, bool doCross) {
    int itr = 0;
    double dl = doCross ? mClosenessTol/10 : -mClosenessTol/10;
    while (itr < mNdtStep) {
        point intp = mDev->intersection(iEdge, ri, rf);
        point r = electron->stepCloseToPoint(intp, dl);
        svec dr = r - intp;

        double d = distance(r, intp);
        if (mDev->intersects(iEdge, r, rf)){
            // did not cross
            if (!doCross && d < mClosenessTol){
                return true;
            }
            electron->doStep();
            ri = r;
            //rf = electron->stepCloseToPoint(intp - dr);

        } else {
            //did cross
            if (doCross && d < mClosenessTol){
                return true;
            }
            rf = r;
            //ri = electron->stepCloseToPoint(intp - dr);
        }
        // FIXME: there might be a better way than this ...
        if (!mDev->intersects(iEdge, ri, rf)) {
            return false;
        }
        itr += 1;
    }
    return false;
}

inline bool Simulator::stepNearEdge2(Particle::ptr electron, point& ri, 
        point& rf, int iEdge, bool doCross) {
    //point rf = intp;
    int itr = 0;
    double dl = doCross ? mClosenessTol : -mClosenessTol;
    point intp = mDev->intersection(iEdge, ri, rf);
    iEdge = doCross ? iEdge : -1;
    while (itr < mNdtStep) {
        intp = electron->stepCloseToPoint(intp, dl);
        if (mDev->intersects(ri, intp) == iEdge) {
            return true;
        }
        itr += 1;
    }
    return false;
}

}}

