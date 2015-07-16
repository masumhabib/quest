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
}

double Particle::timeToReach(const svec& pos) {
    double t = sqrt(pow(pos[0]-r[0],2) + pow(pos[1]-r[1], 2))/speed;
    return t;
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
	//TODO calculation can be made faster if manual calculation is used instead
	// of matrix multiplication
	mat rotationMat;
	// Building rotation matrix
	rotationMat << cos(thti) << -sin(thti) << endr
					<< sin(thti) <<  cos(thti) <<  endr;
	// converting row vector to column vector for matrix multiplication
	col vcol = conv_to< col >::from(this->v);
	col rotv = rotationMat * vcol;
	this->v = conv_to< svec >::from(rotv);
}

}}





