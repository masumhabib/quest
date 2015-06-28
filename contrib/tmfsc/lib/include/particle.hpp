/** Defines the behavior of a praticle
 *  @author Masum Habib<masumhabib.com>
 */

#ifndef TMFSC_PARTICLE_H
#define TMFSC_PARTICLE_H

#include "tmfsc.h"

namespace qmicad { namespace tmfsc {

class Particle {
public:
    virtual const svec& getPos() const { return r; };
    virtual const svec& getVel() const { return v; };
    virtual const svec& getAcc() const { return a; };
    virtual double getTimeStep() const { return dt; };

    virtual void setPotEnergy(double newV) { V = newV; update(); };
    virtual void setMagField(const svec& newB) { B = newB; update(); };
    virtual void setEnergy(double newEn) { En = newEn; update(); };
    virtual void setTimeStep(double dt) { this->dt = dt; update(); };


    virtual const svec& nextPos() = 0;
    virtual void doStep() = 0;
    virtual void reflect(const svec& normal);
    virtual double timeToReach(const svec& pos);

protected:
    virtual void update() = 0;
    Particle(const svec& ri, const svec& vi, double m, double q);

protected:
    double m, q;
    svec r, v;
    svec a = {0,0};
    double speed = 0;

    double V = 0;
    svec B = {0,0,0};
    double En = 0;

    double dt = 0;
    svec nxtr = {0,0};
    svec nxtv = {0,0};

};

}}

#endif



