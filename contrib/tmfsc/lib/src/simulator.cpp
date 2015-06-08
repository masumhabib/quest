/**
 * file: simulator.cpp
 * author: K M Masum Habib
 */

#include "simulator.h"

namespace qmicad { namespace tmfsc {

Simulator::Simulator(Device &dev)
:mDev(dev) {
    mNSteps = 10000;
    mPtsPerCycle = 100;
    mNdtStep = 1000;
    mvF = 1E6;

    mdl = 50;
    mNth = 50;
}

inline tuple<svec, svec, double> Simulator::doStep(const svec &vi, double thi, 
        const svec &ri, double dth, double dt) {
    double thf = thi + dth;
    svec vf = {mvF*cos(thf), mvF*sin(thf)};
    svec rf = ri + (vi + vf)/2*dt;
    return make_tuple(vf, rf, thf);
}

vector<point> Simulator::calcTraj(point ri, double thi, double B, 
            double EF, double V, bool saveTraj) {
    if (fabs(EF-V) < ETOL) {
        throw invalid_argument("EF and V are too close.");
    }

    svec vi = {mvF*cos(thi), mvF*sin(thi)}; // inital velocity
    double wc = mvF*mvF*B/(EF-V); //cyclotron frequency
    double dt = abs(2*pi/wc/mPtsPerCycle); // time step in cyclotron cycle
    double dth = wc*dt; // angle step in cyclotron cycle
    point rf = ri;
    svec vf = vi;
    double thf = thi;
 
    vector<point> r; // trajectory
    if (saveTraj) {
        r.push_back(ri/AA);
    }

    int ii = 0;
    while (ii < mNSteps) {
        cout << "ii " << ii << endl;
        cout << "ri " << ri/nm << endl;
        // get the next step
        tie(vf, rf, thf) = doStep(vi, thi, ri, dth, dt);
        cout << "rf " << rf/nm << endl;

        // check if electron crossed an edge
        int iEdge = mDev.intersects(ri, rf);
        if (iEdge != -1) { 
            // crossed an edge, get the intersection point
            point intp = mDev.intersection(iEdge, ri, rf);

            if (mDev.isReflectEdge(iEdge)) {
                // this is a reflective edge, find the intersection point 
                // and get the time it needs to reach that point
                double dt2 = sqrt(pow(intp[0]-ri[0],2) 
                        + pow(intp[1]-ri[1], 2))/mvF;

                // see if we are close engough to the edge, if not 
                // chenge the time continue until we are inside the device
                // FIXME: the following loop might need optimization.
                while (dt2 > 0) {
                    double dth2 = wc*dt2;
                    tie(vf, rf, thf) = doStep(vi, thi, ri, dth2, dt2);
                    iEdge = mDev.intersects(ri, rf);
                    if (iEdge == -1) {
                        break;
                    }
                    dt2 = dt2 - dt/mNdtStep;
                }

                // reflect electron by changing their velocity direction
                svec a1 = dot(vf, mDev.edgeUnitVect(iEdge))*mDev.edgeUnitVect(iEdge);
                svec a2 = vf - a1;
                vf = a1 - a2;
                thf = atan2(vf[1], vf[0]);
            } else if (mDev.isAbsorbEdge(iEdge)) {
                putElectron(iEdge);
                break;
            }
        }    

        // reset ourselves, ready for the next step
        if (saveTraj) {
            r.push_back(rf/AA);
        }
        ri = rf;
        thi = thf;
        vi = vf;
        ii += 1;
    }

    // last point that we have missed.
    if (saveTraj) {
        r.push_back(rf/AA);
    }
 
    return r;
}

mat Simulator::calcTran(double B, double E, double V, int injCont){
    int nc = mDev.numConts();
    if (injCont < 0 || injCont >= nc){
        throw invalid_argument("Contact number out of bounds");
    }
        
    vector<point> injPts = mDev.createPointsOnCont(injCont, mdl);
    double th0 = mDev.contDirctn(injCont);

    // reset the bins
    resetElectCounts();
    int ne = 0;
    int npts = injPts.size();

    for (int ip = 0; ip < npts; ip += 1) {
        point ri = injPts[ip]*AA;
        vector<double> th(mNth);
        genNormalDist(th, 0, pi/5);

        for (double thi:th){
            if (abs(thi) < pi/2.0-pi/20.0) {
                vector<point> r = calcTraj(ri, th0 + thi + pi/2, B, E, V, 
                        mSaveTraj);
                if (mSaveTraj) {
                    mTraj.push_back(r);
                }
                ne += 1;
            }
        }
        if (mShowProgress == true) {
            int percent = ((ip+1)*100)/npts; 
            cout << percent << "%" << endl;
        }
    }

    mat TE = zeros<mat>(nc,nc);

    for (int ic = 0; ic < nc; ic += 1) {
        int ie = mDev.contToEdgeIndx(ic);
        TE(injCont, ic) = double(mElects[ie])/double(ne);
        TE(ic, injCont) = double(mElects[ie])/double(ne);
    }
    return TE;
}

void Simulator::putElectron(int iEdge, int n){
    mElects[iEdge] += n;
}

void Simulator::resetElectCounts() {
    mElects.resize(mDev.numEdges());
    for (int i = 0; i < mElects.size(); i += 1) {
        mElects[i] = 0;
    }
}



}}

