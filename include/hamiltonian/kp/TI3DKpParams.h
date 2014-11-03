/* 
 * File:   TI3DKpParams.h
 * Author: K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on September 12, 2014, 9:49 AM
 */

#ifndef TI3DKPPARAMS_H
#define	TI3DKPPARAMS_H


#include "maths/constants.h"
#include "maths/arma.hpp"
#include "utils/Printable.hpp"
#include "atoms/AtomicStruct.h"
#include "hamiltonian/hamiltonian.hpp"

namespace qmicad{namespace hamiltonian{

using namespace maths::armadillo;
using namespace maths::constants;

class TI3DKpParams: public cxhamparams{
public:
    TI3DKpParams(const string material="Bi2Se3", const string &prefix = "");
    
    double a() const {return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double A1() const {return mA1; }
    void   A1(double newA1 ){ mA1 = newA1; update(); }
    double A2() const {return mA2; }
    void   A2(double newA2 ){ mA2 = newA2; update(); }
    double B1() const {return mB1; }
    void   B1(double newB1 ){ mB1 = newB1; update(); }
    double B2() const {return mB2; }
    void   B2(double newB2 ){ mB2 = newB2; update(); }
    double C() const {return mC; }
    void   C(double newC ){ mC = newC; update(); }
    double D1() const {return mD1; }
    void   D1(double newD1 ){ mD1 = newD1; update(); }
    double D2() const {return mD2; }
    void   D2(double newD2 ){ mD2 = newD2; update(); }
    double M() const {return mM; }
    void   M(double newM ){ mM = newM; update(); }

    virtual string toString() const;

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    
private:
    //!< Default parameters.
    void setDefaultParams(){
        ma     = 5.0;
        mdtol  = 1E-3; 
        mortho = true;
        mpt.add(0, "D", 4, 4);
        
    };  
    
    void setBi2Se3Params(){
        // Bi2Se3 parameters from: Nat Phys 9, 438 (2009). 
        mA1    = 2.2;       // eV*A
        mA2    = 4.1;       // eV*A
        mB1    = 10;        // eV*A^2
        mB2    = 56.6;      // eV*A^2
        mC     = -0.0068;   // eV
        mD1    = 1.3;       // eV*A^2
        mD2    = 19.6;      // eV*A^2
        mM     = 0.28;      // eV

    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update();

private:
    //!< k.p parameters
    double ma;            // discretization length in y direction
    //!< k.p parameters, see Rev. Mod. Phys. 83, 1057 (2011) or Nat Phys 9, 438 (2009).
    double mA1;           // A1
    double mA2;           // A2
    double mB1;           // B1
    double mB2;           // B2
    double mC;            // C 
    double mD1;           // D1
    double mD2;           // D2
    double mM;            // M 
    
    //!< tight binding parameters
    cxmat mI;
    cxmat mU;            // unitary transformation matrix
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    cxmat mt01z;
    cxmat mt10z;    
    
};

}}
#endif	/* TI3DKPPARAMS_H */

