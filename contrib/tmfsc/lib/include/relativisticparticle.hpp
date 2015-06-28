/**
 * @author K M Masum Habib <masumhabib.com>
 */

#ifndef TMFSC_RELATIVISTICPARTICLE_HPP
#define TMFSC_RELATIVISTICPARTICLE_HPP

#include "particle.hpp"
#include "maths/constants.h"

namespace qmicad { namespace tmfsc {
using maths::constants::e;

class RelativisticParticle : public Particle {
public:
    RelativisticParticle(const svec& ri, const svec& vi, double En, 
            double V = 0, double B = 0, double q = e);
    virtual const svec& nextPos();
    virtual void doStep();
    virtual void reflect(const svec& normal);
    virtual double timeToReach(const svec& pos);

    void setFermiVel(double newVF) { vF = newVF; update(); };

protected:
    virtual void update();

protected:
    double vF;
    double th = 0;
    double nxtth = 0, dth = 0;
    double wc;
};

}}


#endif





