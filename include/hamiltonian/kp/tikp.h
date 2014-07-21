/* 
 * File:   ti.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef TIKP_H
#define	TIKP_H

#include <armadillo>

#include "maths/constants.h"
#include "utils/Printable.hpp"
#include "atoms/AtomicStruct.h"
#include "hamiltonian/hamiltonian.hpp"

namespace qmicad{
namespace hamiltonian{

using boost::static_pointer_cast;
using namespace maths::armadillo;
using namespace maths::constants;

/**
 * TI surface k.p Hamiltonian parameters with two spins per site.
 * In two spin basis [z_up, z_dn], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud
 *  H_du  Hdd
 * A quantity, sigma_z*(kx^2 + ky^2) is added to the hamiltonian to avoid 
 * the fermion doubling problem. 
 * See: Phys. Rev. B 86, 085131 (2012) and the references therein.
 */
class TISurfKpParams: public cxhamparams{
public:    
    TISurfKpParams(const string &prefix = "");
        
    double a(){return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double C(){return mC; }
    void   C(double newC ){ mC = newC; update(); }
    double A2(){return mA2; }
    void   A2(double newA2 ){ mA2 = newA2; update(); }
    double K(){return mK; }
    void   K(double newK ){ mK = newK; update(); }
    
    virtual string toString() const;

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const;
    
private:
    //!< Default parameters.
    void setDefaultParams(){
        ma     = 2.0;
        mK     = 1.0;  
        mC     = 0.0; 
        mA2    = 3.33; 
        mdtol  = 1E-3; 
        mortho = true;
    };    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update();
    
private:
    //!< k.p parameters
    double ma;           // discretization length in x direction
    double mK;           // parameter to avoid fermion doubling
    double mC;           // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mA2;          // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    
    //!< tight binding parameters
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    
};

}
}
#endif	/* TIKP_H */

