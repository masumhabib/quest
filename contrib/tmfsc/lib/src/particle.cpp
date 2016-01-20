/**
 * @author K M Masum Habib <masumhabib.com
 * @co-author Mirza Elahi <mirza.monzur@gmail.com
 */

#include "particle.hpp"

namespace qmicad { namespace tmfsc {

Particle::Particle(const svec& ri, const svec& vi, double m, double q) 
        :r(ri), v(vi), m(m), q(q) 
{
    speed = sqrt(v[0]*v[0] + v[1]*v[1]);
    speed2 = speed*speed;
}

double Particle::timeToReach(const svec& pos) {
    double dx = pos[0]-r[0], dy = pos[1]-r[1];
    return sqrt(dx*dx + dy*dy)/speed;
}

const svec& Particle::stepCloseToPoint(const svec& pos, double distanceTol) {
    double dtOld = getTimeStep();
    double timeTol = distanceTol/speed;
    // get the time we need to reach the point
    setTimeStep(timeToReach(pos) + timeTol);
    // got to that point
    nextPos();
    // reset time step
    setTimeStep(dtOld);

    return nxtr;
}

void Particle::reflect(const svec& normal) {
    svec vparallel =  dot(v, normal)*normal;
    svec vparpendicular = v - vparallel;
    v = vparpendicular - vparallel;
}

void Particle::rotateVel(double thti){
	double vx = this->v[0], vy = this->v[1];
	this->v[0] = vx*cos(thti) - vy*sin(thti);
	this->v[1] = vx*sin(thti) + vy*cos(thti);
}

}}





