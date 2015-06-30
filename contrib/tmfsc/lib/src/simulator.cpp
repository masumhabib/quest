/**
 * file: simulator.cpp
 * author: K M Masum Habib
 */

#include "simulator.h"

namespace qmicad { namespace tmfsc {

Simulator::Simulator(Device &dev)
:mDev(dev) {
    std::cout << std::endl << " Simu " << std::endl;
}

Trajectory Simulator::calcTraj(Particle& particle, bool saveTraj) {
    point ri = particle.getPos();
    svec vi = particle.getVel();
    svec rf = ri;

    Trajectory r;
    if (saveTraj) {
        r.push_back(ri);
    }

    int ii = 0;
    while (ii < mMaxStepsPerTraj) {
        // get the next position
        rf = particle.nextPos();

        // check if electron crossed an edge
        int iEdge = mDev.intersects(ri, rf);
        if (iEdge == -1) { // no crossing, continue
            particle.doStep();
        } else { // crossing, check wchich boundary is crossed
            // crossed an edge, get the intersection point
            point intp = mDev.intersection(iEdge, ri, rf);
            if (mDev.isReflectEdge(iEdge)) {
                // this is a reflective edge, find the intersection point 
                // and get the time it needs to reach that point
                double dt = particle.getTimeStep();
                double dt2 = particle.timeToReach(intp)*0.999;

                // see if we are close engough to the edge, if not 
                // chenge the time continue until we are inside the device
                // FIXME: the following loop might need optimization.
                int idtStep = 0;
                while (dt2 > 0 && idtStep < mNdtStep) {
                    particle.setTimeStep(dt2);
                    rf = particle.nextPos();
                    int iEdge2 = mDev.intersects(ri, rf);
                    if (iEdge2 == -1) {
                        particle.doStep();
                        break;
                    }
                    dt2 = dt2 - dt/mNdtStep;
                    idtStep += 1;
                }

                // reflect electron by changing their velocity direction
                particle.setTimeStep(dt);
                particle.reflect(mDev.edgeNormVect(iEdge));
            } else if (mDev.isAbsorbEdge(iEdge)) {
                collectElectron(mDev.edgeToContIndx(iEdge));
                break;
            }
        }    

        // reset ourselves, ready for the next step
        if (saveTraj) {
            r.push_back(rf);
        }
        ri = rf;
        ii += 1;
    }

    // last point that we have missed.
    if (saveTraj) {
        r.push_back(rf);
    }
 
    return r;
}

Trajectory Simulator::calcTraj(point ri, double thi, double B, 
            double EF, double V, bool saveTraj) {
    if (fabs(EF-V) < ETOL) {
        throw invalid_argument("EF and V are too close.");
    }

    svec vi = {mvF*cos(thi), mvF*sin(thi)}; // inital velocity
    if (isAutoDt) {
        double wc = mvF*nm*mvF*nm*B/(EF-V); //cyclotron frequency
        dt = abs(2*pi/wc/mPtsPerCycle); // time step in cyclotron cycle
        if (dt > 1E-3) {
            std::cout << "W:  Too big time step" << std::endl;
        }
    }

    shared_ptr<Particle> particle;
    if (particleType == ParticleType::DiracCyclotron) {
        particle = make_shared<DiracCyclotron>(ri, vi, EF, V, B);
    } else {
        particle = make_shared<DiracElectron>(ri, vi, EF, V, B);
    }
    particle->setTimeStep(dt);
 
    return calcTraj(*particle, saveTraj);
}

tuple<mat, TrajectoryVect> Simulator::calcTran(double B, double E, double V, 
        int injCont, bool saveTraj){
    int nc = mDev.numConts();
    if (injCont < 0 || injCont >= nc){
        throw invalid_argument("Contact number out of bounds");
    }
        
    vector<point> injPts = mDev.createPointsOnCont(injCont, mdl);
    double th0 = mDev.contDirctn(injCont) + pi/2;
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
                Trajectory r = calcTraj(ri, th0 + thi, B, E, V, 
                        saveTraj);
                if (saveTraj) {
                    trajs.push_back(r);
                }
            }
        }
    }

    mat TE = zeros<mat>(nc,nc);

    for (int ic = 0; ic < nc; ic += 1) {
        TE(injCont, ic) = mElectBins[ic]/mnElects;
        TE(ic, injCont) = mElectBins[ic]/mnElects;
    }
    return make_tuple(TE, trajs);
}

void Simulator::collectElectron(int iCont, double n){
    if (iCont < mElectBins.size()) {
        mElectBins[iCont] += n;
        mnElects += n;
    }
}

void Simulator::resetElectBins() {
    mElectBins.resize(mDev.numConts());
    for (auto &bin : mElectBins) {
        bin = 0;
    }
    mnElects = 0;
}



}}

