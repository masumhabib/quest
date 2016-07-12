/**
 * @author K M Masum Habib <masumhabib.com>
 */

#ifndef TMFSC_DIRACCYCLOTRON_HPP
#define TMFSC_DIRACCYCLOTRON_HPP

#include "particle.hpp"
#include "maths/constants.h"
#include "maths/arma.hpp"
#include <memory>

namespace qmicad { namespace tmfsc {
using maths::constants::e;
using std::make_shared;

class DiracCyclotron : public Particle {
public:
    DiracCyclotron(const svec& ri, const svec& vi, double En, 
            double V = 0, double Bz = 0);
    virtual const svec& nextPos();
    virtual void doStep();
    virtual void reflect(const svec& normal);
    virtual void refract(const svec& normal, double V1, double V2);
    virtual void rotateVel(double thti);
    virtual ptr clone();

protected:
    virtual void update();

protected:
    double th = 0;
    double nxtth = 0, dth = 0;
    double wc;
    
};

}}


#endif





