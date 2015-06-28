/**
 * @author K M Masum Habib <masumhabib.com>
 *
 */

#include "relativisticparticle.hpp"

namespace qmicad { namespace tmfsc {
RelativisticParticle::RelativisticParticle(const svec& ri, const svec& vi, 
        double En, double V, double B, double q) : Particle(ri, vi, 0, q)
{
    this->En = En;
    this->V = V;
    this->B = B;
    th = atan2(v[1], v[0]);
    vF = sqrt(v[0]*v[0] + v[1]*v[1]);
}

const svec& RelativisticParticle::nextPos() {
    nxtth = th + dth;
    nxtv = svec({vF*cos(nxtth), vF*sin(nxtth)});
    nxtr = r + (v + nxtv)/2*dt;

    return nxtr;
}

void RelativisticParticle::doStep() {
    v = nxtv;
    r = nxtr;
    th = nxtth;
}

void RelativisticParticle::reflect(const svec& normal) {
    svec vparallel =  dot(v, normal)*normal;
    svec vparpendicular = v - vparallel;
    v = vparpendicular - vparallel;

    th = atan2(v[1], v[0]);
}

double RelativisticParticle::timeToReach(const svec& pos) {
    double t = sqrt(pow(pos[0]-r[0],2) + pow(pos[1]-r[1], 2))/vF;
    return t;
}

void RelativisticParticle::update() {
    wc = vF*nm*vF*nm*B/(En-V); //cyclotron frequency
    dth = wc*dt; // angle step in cyclotron cycle
}



}}




