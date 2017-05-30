/** Defines the behavior of a praticle
 *  @author Masum Habib<masumhabib.com>
 */

#ifndef TMFSC_PARTICLE_H
#define TMFSC_PARTICLE_H
#include "tmfsc.h"
#include <memory>
#include "maths/constants.h"

namespace quest { namespace tmfsc {

using maths::armadillo::mat;
using maths::armadillo::col;
using arma::conv_to;
using arma::endr;
using std::shared_ptr;
using maths::constants::pi;

class Particle {
public:
    typedef shared_ptr<Particle> ptr;

    virtual const svec& getPos() const { return r; };
    virtual const svec& getVel() const { return v; };
    virtual const svec& getAcc() const { return a; };
    virtual double getTimeStep() const { return dt; };
    virtual void setTimeStep(double dt) { this->dt = dt; update(); };
    virtual double getPot() const { return V; };
    virtual void setPot(double newV) { V = newV; update(); };
    virtual void setMagField(const svec& newB) { B = newB; update(); };
    virtual double getEnergy() const { return En; };
    virtual void setEnergy(double newEn) { En = newEn; update(); };
    
    virtual void setOccupation(double newOcc) { occupation = newOcc; };
    virtual double getOccupation() const { return occupation; };

    virtual const svec& nextPos() = 0;
    virtual void doStep() = 0;
    virtual ptr clone() = 0;
    virtual void reflect(const svec& normal);
    virtual void refract(const svec& normal, double V1, double V2);
    virtual void rotateVel(double thti); //!< Rotate the Particle by thti
    virtual double timeToReach(const svec& pos);
    virtual const svec& stepCloseToPoint(const svec& pos, 
            double distanceTol = 0.0);

    virtual void flipCharge();
    //operators
    friend bool operator< (const Particle& lhs, const Particle& rhs){ 
        return lhs.occupation < rhs.occupation; 
    }
    friend bool operator> (const Particle& lhs, const Particle& rhs) {return rhs < lhs;}
    friend bool operator<=(const Particle& lhs, const Particle& rhs) {return !(lhs > rhs);}
    friend bool operator>=(const Particle& lhs, const Particle& rhs) {return !(lhs < rhs);}

protected:
    virtual void update() = 0;
    Particle(const svec& ri, const svec& vi, double m, double q);

protected:
    double m, q;
    svec r, v;
    svec a = {0,0};
    double speed = 0;
    double speed2 = 0; // speed^2
    double occupation = 1.0;

    double V = 0;
    svec B = {0,0,0};
    double En = 0;

    double dt = 0;
    svec nxtr = {0,0};
    svec nxtv = {0,0};

};

// Comparator for smart ptr of Particle
struct ParticleComparator {
    bool operator() (Particle::ptr lhs, Particle::ptr rhs) {
        return *lhs < *rhs;
    }
};

}}

#endif



