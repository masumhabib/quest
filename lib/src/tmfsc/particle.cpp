/**
 * @author K M Masum Habib <masumhabib.com
 * @co-author Mirza Elahi <mirza.monzur@gmail.com
 */

#include "tmfsc/particle.hpp"

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

void Particle::refract(const svec& normal, double V1, double V2) {
	double thti, ratio, costhetai, sin2thetar;
	svec v_refrac, v_unit;
	double vMag = sqrt( v[0]*v[0] + v[1]*v[1] );
	v_unit = v/vMag;
	svec NormVec = normal;
	// incident angle
	thti = acos( dot( NormVec, v_unit ) );
	// correction for incident angle for normal direction not in traditional
	// direction
	if( thti > pi/2 ){
		thti = pi - thti;
	}else{
		NormVec = -NormVec;
	}
	ratio = V1/V2;
	costhetai = cos(thti);
	sin2thetar = ratio * ratio * ( 1- costhetai*costhetai );
	v_refrac = ratio * v_unit + (ratio*costhetai - sqrt(1-sin2thetar)) * NormVec;
	// reference
	// http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
	v = v_refrac * vMag;
}

void Particle::rotateVel(double thti){
	double vx = this->v[0], vy = this->v[1];
	this->v[0] = vx*cos(thti) - vy*sin(thti);
	this->v[1] = vx*sin(thti) + vy*cos(thti);
}

}}





