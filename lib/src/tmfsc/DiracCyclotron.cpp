/**
 * @author K M Masum Habib <masumhabib.com>
 * @co-author Mirza Elahi <mirza.monzur@gmail.com>
 */

#include "tmfsc/DiracCyclotron.hpp"

namespace quest { namespace tmfsc {
DiracCyclotron::DiracCyclotron(const svec& ri, const svec& vi, 
        double En, double V, double Bz) : Particle(ri, vi, 0, -e)
{
    this->En = En;
    this->V = V;
    B = {0, 0, Bz};
    th = atan2(v[1], v[0]);
    update();
}

const svec& DiracCyclotron::nextPos() {
    nxtth = th + dth;
    nxtv = svec({speed*cos(nxtth), speed*sin(nxtth)});
    nxtr = r + (v + nxtv)/2*dt;
    return nxtr;
}

void DiracCyclotron::doStep() {
    v = nxtv;
    r = nxtr;
    th = nxtth;
}

void DiracCyclotron::reflect(const svec& normal) {
    Particle::reflect(normal);
    th = atan2(v[1], v[0]);

}

void DiracCyclotron::refract(const svec& normal, double V1, double V2){
	Particle::refract(normal, V1, V2);
	th = atan2(v[1], v[0]);
}

void DiracCyclotron::rotateVel(double thti) {
    Particle::rotateVel(thti);

    th = atan2(v[1], v[0]);
}


DiracCyclotron::ptr DiracCyclotron::clone() {
    auto newElect = make_shared<DiracCyclotron>(*this);
    return newElect;
}

void DiracCyclotron::update() {
    double Bz = B[2];
    int sign_q = ( q >= 0 ? 1 : -1 ); // sign of q
    wc = - sign_q * speed2*nm2*Bz/(En-V); // cyclotron frequency
    dth = wc*dt; // angle step in cyclotron cycle
}

void DiracCyclotron::flipCharge(){
    Particle::flipCharge();
    update();
}

}}




