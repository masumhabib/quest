/**
 * @author K M Masum Habib <masumhabib.com>
 */

#ifndef TMFSC_DIRACELECTRON_HPP
#define TMFSC_DIRACELECTRON_HPP

#include "particle.hpp"
#include "maths/constants.h"
#include "maths/arma.hpp"
#include <memory>

namespace qmicad { namespace tmfsc {
using maths::constants::e;
using std::make_shared;

class DiracElectron : public Particle {
public:
    DiracElectron(const svec& ri, const svec& vi, double En, 
            double V = 0, double Bz = 0);
    virtual const svec& nextPos();
    virtual void doStep();
    virtual void reflect(const svec& normal);
    virtual void refract(const svec& normal, double n1, double n2);
    virtual ptr clone();
    virtual void flipCharge();

protected:
    virtual void update();
    void updateAcceleration();
    inline void updateLorentzForce();

protected:
    svec F;
};

}}


#endif





