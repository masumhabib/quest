/* 
 * File:   tikp4.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#ifndef TIKP4_H
#define	TIKP4_H

#include "maths/constants.h"
#include "maths/arma.hpp"
#include "utils/Printable.hpp"
#include "atoms/AtomicStruct.h"
#include "hamiltonian/hamiltonian.hpp"

namespace quest{
namespace hamiltonian{

using namespace maths::armadillo;
using namespace maths::constants;

/**
 * TI surface k.p Hamiltonian parameters in four spins basis set.
 * In four spin basis [z_up, z_dn, z_up, z_dn], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud  0    0
 *  H_du  Hdd  0    0
 *  0     0    Huu  Hud
 *  0     0    Hdu  Hdd
 * A quantity, sigma_z*(kx^2 + ky^2) is added to the upper 2x2 and subtracted 
 * from the low 2x2 to avoid fermion doubling and to avoid all the side
 * effects that can be introduced by the term sigma_z*(kx^2 + ky^2).
 */
class TISurfKpParams4: public cxhamparams{
    
public:
    TISurfKpParams4(const string &prefix = ""); 
    
    double a() const {return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double C() const {return mC; }
    void   C(double newC ){ mC = newC; update(); }
    double A2() const {return mA2; }
    void   A2(double newA2 ){ mA2 = newA2; update(); }
    double K() const {return mK; }
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
        mpt.add(0, "D", 4, 4);
    };    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update();
    
private:
    //!< k.p parameters
    double ma;            // discretization length in y direction
    double mK;            // parameter to avoid fermion doubling
    double mC;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mA2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    
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
#endif	/* TIKP4_H */

