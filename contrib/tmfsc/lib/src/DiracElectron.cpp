/**
 * @author K M Masum Habib <masumhabib.com>
 *
 */

#include "DiracElectron.hpp"

namespace qmicad { namespace tmfsc {
DiracElectron::DiracElectron(const svec& ri, const svec& vi, 
        double En, double V, double Bz) : Particle(ri, vi, 0, -e)
{
    this->En = En;
    this->V = V;
    B = {0,0,Bz};
    update();
}

const svec& DiracElectron::nextPos() {
    nxtv = v + a*dt;
    nxtr = r + v*dt;
    return nxtr;
}

void DiracElectron::doStep() {
    v = nxtv;
    r = nxtr;
    update();
}

void DiracElectron::reflect(const svec& normal) {
    Particle::reflect(normal);
    update();
}

void DiracElectron::refract(const svec& normal, double n1, double n2){
	Particle::refract(normal, n1, n2);
	update();
}

DiracElectron::ptr DiracElectron::clone() {
    return make_shared<DiracElectron>(*this); 
}

void DiracElectron::update() {
    double Bz = B[2];
    F = {q*v[1]*Bz, -q*v[0]*Bz};
    //speed = sqrt(v[0]*v[0] + v[1]*v[1]);
    m = e*(En - V)/(speed*speed);
    a = F/m*nm*nm;
}

}}




