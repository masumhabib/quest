/**
 * @author K M Masum Habib <masumhabib.com
 */

#include "particle.hpp"

namespace qmicad { namespace tmfsc {

Particle::Particle(const svec& ri, const svec& vi, double m, double q) 
        :r(ri), v(vi), m(m), q(q) 
{
    speed = sqrt(v[0]*v[0] + v[1]*v[1]);
}

double Particle::timeToReach(const svec& pos) {
    double t = sqrt(pow(pos[0]-r[0],2) + pow(pos[1]-r[1], 2))/speed;
    return t;
}

void Particle::reflect(const svec& normal) {
    svec vparallel =  dot(v, normal)*normal;
    svec vparpendicular = v - vparallel;
    v = vparpendicular - vparallel;
}


}}





